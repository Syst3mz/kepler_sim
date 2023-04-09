use macroquad::input::KeyCode::{Enter, LeftControl, RightControl};
use macroquad::main;
use macroquad::miniquad::KeyCode::RightShift;
use macroquad::prelude::*;
use macroquad::prelude::KeyCode::LeftShift;
use num_traits::FloatConst;

struct State {
    orbiters: Vec<Satellite>,
    orbiters_last_pos: Vec<DVec2>,
    global_time: f64,
    global_time_scale: f64
}

impl State {
    fn update(&mut self) {
        // update time
        self.global_time += (get_frame_time() as f64) * self.global_time_scale;
    }
    fn render(&self) {
        draw_text(&format!("Timescale: {}", self.global_time_scale), 0.0, 20.0, 15.0, WHITE);
        for orbiter in &self.orbiters {
            let position = Satellite::compute_point_on_orbit_2d(
                orbiter.apoapsis,
                orbiter.periapsis,
                orbiter.global_time_to_true_anomaly(self.global_time));
            let text = format!("dist: {:.2}", DVec2::distance(position, DVec2 {x: 0.0, y: 0.0}));

            draw_orbit_lines(orbiter.offset.x,
                             orbiter.offset.y,
                             orbiter.apoapsis,
                             orbiter.periapsis,
                             360, 3.0, GRAY
            );

            draw_circle(position.x as f32 + orbiter.offset.x,
                        position.y as f32 + orbiter.offset.y,
                        20.0 * orbiter.scale, WHITE);
            draw_text(&text,
                      position.x as f32 + orbiter.offset.x,
                      position.y as f32 + orbiter.offset.y,
                      12.0, RED);
            // 0 0, sm
        }
    }
}

fn draw_orbit_lines(center_x: f32, center_y: f32, apoapsis: f64, periapsis: f64, segments: u32, thickness: f32, color: Color) {
    let angle_step = 2.0 * (f64::PI() / segments as f64);
    for segment_number in 0..segments {
        let first_segment = Satellite::compute_point_on_orbit_2d(apoapsis, periapsis, segment_number as f64 * angle_step);
        let second_segment = Satellite::compute_point_on_orbit_2d(apoapsis, periapsis, (segment_number + 1) as f64 * angle_step);

        draw_line(first_segment.x as f32 + center_x,
                  first_segment.y as f32 + center_y,
                  second_segment.x as f32 + center_x,
                  second_segment.y as f32 + center_y, thickness, color)
    }
}

#[main("Orbit Simulation")]
async fn main() {
    let mut state = setup();
    println!("working in a {}x{} window.", screen_width(), screen_height());
    loop {
        clear_background(BLACK);
        draw_circle(half_screen_width(), half_screen_height(), 45.0, YELLOW);
        //draw_ellipse_lines(half_screen_width(), half_screen_height(), 100.0, 50.0, 24, 3.0, GRAY);

        handle_input(&mut state);

        state.update();
        state.render();
        if get_frame_time() > 1.0 / 60.0 {
            println!("Fame took too long: {}", get_frame_time())
        }

        next_frame().await
    }
}

fn handle_input(state: &mut State) {

    let mut speed_delta_factor = 1.0;
    if is_key_down(LeftShift) || is_key_down(RightShift) {
        speed_delta_factor *= 10.0;
    }
    if is_key_down(KeyCode::Z) {
        speed_delta_factor *= 100.0;
    }

    if is_key_pressed(KeyCode::L) {
        state.global_time_scale += 1.0 * speed_delta_factor;
    }
    else if is_key_pressed(KeyCode::J) {
        state.global_time_scale -= 1.0 * speed_delta_factor;
    } else if is_key_pressed(KeyCode::K) {
        state.global_time_scale = 10.0
    }
}

fn half_screen_width() -> f32 {
    screen_width() / 2.0
}

fn half_screen_height() -> f32 {
    screen_height() / 2.0
}

fn setup() -> State {
    State {
        orbiters: vec![
            Satellite {
                apoapsis: 200.0,
                periapsis: 200.0,
                scale: 1.0,
                offset: Vec2 { x: half_screen_width(), y: half_screen_height() }
            },

            Satellite {
                apoapsis: 250.0,
                periapsis: 150.0,
                scale: 0.5,
                offset: Vec2 { x: half_screen_width(), y: half_screen_height() },
            },

            Satellite {
                apoapsis: 150.0,
                periapsis: 100.0,
                scale: 0.25,
                offset: Vec2 { x: half_screen_width(), y: half_screen_height() },
            },
            Satellite {
                apoapsis: 300.0,
                periapsis: 50.0,
                scale: 0.25,
                offset: Vec2 { x: half_screen_width(), y: half_screen_height() },
            }
        ],
        orbiters_last_pos: vec![DVec2::default(), DVec2::default(), DVec2::default()],
        global_time: 0.0,
        global_time_scale: 10.0
    }
}

// M = E - esin(E)

const EPSILON: f64 = 0.000001;
fn newtons_method<F, G>(x0: f64, fx: F, dx: G) -> f64
    where F: Fn(f64) -> f64, G: Fn(f64) -> f64 {
    let mut root = x0;

    loop {
        let next_root = root - (fx(root) / dx(root));
        let delta = f64::abs(next_root - root);
        root = next_root;

        if delta <= EPSILON {
            break;
        }
    }

    root
}

const CENTRAL_BODY_MASS:f64 = 10.0;
trait Orbiter {

    fn keplers_equation(mean_anomaly: f64, eccentric_anomaly: f64, eccentricity: f64) -> f64{
        return eccentric_anomaly - (eccentricity * f64::sin(eccentric_anomaly)) - mean_anomaly
    }

    fn keplers_equation_derivation(eccentric_anomaly: f64, eccentricity: f64) -> f64 {
        return 1.0 -(eccentricity * f64::cos(eccentric_anomaly))
    }

    fn solve_keplers(mean_anomaly: f64, eccentricity: f64) -> f64 {
        let eccentric_anomaly = if eccentricity > 0.8 { f64::PI() } else {mean_anomaly};
        newtons_method(
            eccentric_anomaly,
            |eccentric_anomaly| Self::keplers_equation(mean_anomaly, eccentric_anomaly, eccentricity),
            |eccentric_anomaly| Self::keplers_equation_derivation(eccentric_anomaly, eccentricity)
        )
    }

    fn compute_orbital_period(semi_major_axis: f64) -> f64 {
        2.0 * f64::PI() * f64::sqrt(f64::powi(semi_major_axis, 3) / CENTRAL_BODY_MASS)
    }

    fn compute_semi_major_axis(apoapsis: f64, periapsis:f64) -> f64 {
        (apoapsis + periapsis) / 2.0
    }
    fn compute_semi_minor_axis(apoapsis: f64, periapsis:f64) -> f64 {
        f64::sqrt(apoapsis * periapsis)
    }

    fn compute_point_on_orbit_2d(apoapsis: f64, periapsis: f64, true_anomaly:f64) -> DVec2 {
        let semi_major_axis = Self::compute_semi_major_axis(apoapsis, periapsis);
        let semi_minor_axis = Self::compute_semi_minor_axis(apoapsis, periapsis);

        let mean_anomaly = true_anomaly * f64::PI() * 2.0;
        let linear_eccentricity = semi_major_axis - periapsis;
        let eccentricity = linear_eccentricity / semi_major_axis;

        let eccentric_anomaly = Self::solve_keplers(mean_anomaly, eccentricity);

        let x = semi_major_axis * (f64::cos(eccentric_anomaly) - eccentricity);
        let y = semi_minor_axis * f64::sin(eccentric_anomaly);

        dvec2(x, y)
    }
}

struct Satellite {
    apoapsis: f64,
    periapsis: f64,
    scale: f32,
    offset: Vec2
}
impl Orbiter for Satellite {}

impl Satellite {
    fn global_time_to_true_anomaly(&self, global_time: f64) -> f64 {
        let orbital_period = Self::compute_semi_major_axis(self.apoapsis, self.periapsis);
        global_time / orbital_period
    }
}


#[cfg(test)]
pub mod test {
    use crate::{EPSILON, newtons_method, Orbiter, Satellite};

    #[test]
    fn check_five_minus_x() {
        let ans = newtons_method(
            0.1,
            |x| 5.0-x,
            |_x| -1.0
        );
        assert!(f64::abs(5.0 - ans) <= EPSILON)
    }

    #[test]
    fn check_x_squared_over_3_minus_3() {
        let ans = newtons_method(
            0.1,
            |x| (f64::powi(x, 2) / 3.0) - 3.0,
            |x| (2.0 * x) / 3.0
        );
        assert!(f64::abs(3.0 - ans) <= EPSILON)
    }
    #[test]
    fn check_eccentricity() {
        let ap = 100.0;
        let per = 2.0;

        let sma = Satellite::compute_semi_major_axis(ap, per);
        let sna = Satellite::compute_semi_minor_axis(ap ,per);
        let lin_ecc = sna - per;
        let ecc = lin_ecc / sma;
        println!("ap: {}, per {}, sma: {}, sna: {}, lin_ecc: {}, ecc: {}", ap, per, sma, sna, lin_ecc, ecc);
        assert!(false)
    }
}