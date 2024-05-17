use macroquad::prelude::*;

const TARGET_FPS: u32 = 60;
const TARGET_DT: f32 = 1. / TARGET_FPS as f32;

#[derive(Clone, Copy)]
/// Positon is not relative to the screen size
/// It is in range [-1, 1]
/// where (0, 0) is the center of the screen
/// and 1 in any direction is half of smallest screen dimension
struct Particle {
    pos: Vec2,
    prev_pos: Vec2,
    color: Color,
}

impl Particle {
    /// The radius of the particle
    /// With a radius of 1, the particle will be as big as the screen
    const RADIUS: f32 = 0.01;

    const GRAVITY: f32 = 0.5;
    const GRAVITY_VEC: Vec2 = Vec2::new(0., Self::GRAVITY);

    const CONSTRAINT_CENTER: Vec2 = Vec2::ZERO;

    fn random_color() -> Color {
        Color::new(
            rand::gen_range(0., 1.),
            rand::gen_range(0., 1.),
            rand::gen_range(0., 1.),
            1.,
        )
    }
    fn new(position: Vec2, initial_velocity: Vec2, color: Option<Color>) -> Self {
        Particle {
            pos: position,
            prev_pos: position - initial_velocity * TARGET_DT,
            color: color.unwrap_or_else(|| Self::random_color()),
        }
    }

    fn update_verlet(&mut self, dt: f32) {
        let new_pos = 2. * self.pos - self.prev_pos + Self::GRAVITY_VEC * dt * dt;
        self.prev_pos = self.pos;
        self.pos = new_pos;

        self.constraint();
    }

    fn draw(&self) {
        let dim = Simulation::dimension();

        let screen_x = screen_width() * 0.5 + self.pos.x * dim;
        let screen_y = screen_height() * 0.5 + self.pos.y * dim;

        draw_circle(
            screen_x,
            screen_y,
            Self::RADIUS * dim,
            self.color,
        );
    }

    fn constraint(&mut self) {
        let delta = self.pos - Self::CONSTRAINT_CENTER;
        let max_dist = Simulation::CONSTRAINT_RADIUS - Self::RADIUS;
        let dist = delta.length();

        if dist > max_dist{
            // self.pos = Self::CONSTRAINT_CENTER + delta.normalize() * (Simulation::CONSTRAINT_RADIUS - Self::RADIUS);
            self.pos -= delta.normalize() * (dist - max_dist);
        }
    }

    fn distance(&self, other: &Particle) -> f32 {
        (self.pos - other.pos).length()
    }
}

struct Simulation {
    particles: Vec<Particle>,
}

impl Simulation {
    /// 90% of the whole area
    const CONSTRAINT_RADIUS: f32 = 0.9 * 0.5;

    fn dimension() -> f32 {
        screen_width().min(screen_height())
    }

    fn new() -> Self {
        Simulation {
            particles: vec![Particle::new(Vec2::new(0., 0.), Vec2::new(0., 0.), None)],
        }
    }

    fn handle_collision(&mut self, a_id: usize, b_id: usize) {
        let a = &self.particles[a_id];
        let b = &self.particles[b_id];

        let dist_a_b = a.distance(b);
        let max_dist = 2. * Particle::RADIUS;
        let collision = dist_a_b < max_dist;

        if a_id == b_id || !collision {
            return;
        }

        let delta_a_b = (b.pos - a.pos).normalize() * (max_dist - dist_a_b) * 0.5;

        let a = &mut self.particles[a_id];
        a.pos -= delta_a_b;

        let b = &mut self.particles[b_id];
        b.pos += delta_a_b;
    }

    fn handle_colliosions(&mut self) {
        for a_id in 0..self.particles.len() {
            for b_id in 0..self.particles.len() {
                self.handle_collision(a_id, b_id);
            }
        }
    }

    fn update(&mut self, dt: f32) {
        for particle in self.particles.iter_mut() {
            particle.update_verlet(dt);
        }
        self.handle_colliosions();
    }

    fn draw_constraint(&self) {
        draw_poly_lines(
            screen_width() * 0.5,
            screen_height() * 0.5,
            64,
            Self::CONSTRAINT_RADIUS * Self::dimension(),
            0.,
            2.,
            WHITE,
        );
    }

    fn draw(&self) {
        for particle in self.particles.iter() {
            particle.draw();
        }
        self.draw_constraint();
    }

    fn spawn_particle(&mut self) {
        let pos = Vec2::new(-0.3, 0.);
        let vel = Vec2::new(0.5, 0.);

        let particle = Particle::new(pos, vel, None);
        self.particles.push(particle);
    }

    fn spawn_at_mouse(&mut self) {
        let mouse = mouse_position();
        let dimension = Self::dimension();

        let x = (mouse.0 - screen_width() * 0.5)/ dimension;
        let y = (mouse.1 - screen_height() * 0.5) / dimension;
        let pos = Vec2::new(x, y);
        let vel = Vec2::ZERO;

        let particle = Particle::new(pos, vel, None);
        self.particles.push(particle);
    }
}

struct GameState {
    simulation: Simulation,
}

impl GameState {

    const FONT_SIZE: f32 = 30.;

    fn new() -> Self {
        GameState {
            simulation: Simulation::new(),
        }
    }

    fn handle_inputs(&mut self) {
        if is_key_pressed(KeyCode::S) {
            self.simulation.spawn_particle();
        }
        if is_key_down(KeyCode::A) {
            self.simulation.spawn_particle();
        }
        if is_mouse_button_pressed(MouseButton::Left) {
            self.simulation.spawn_at_mouse();
        }
    }

    fn display_stats(&self) {
        let particles = self.simulation.particles.len();
        let text = format!("Particles: {}", particles);

        draw_text(&text, 10., 20., Self::FONT_SIZE, WHITE);

        let fps = get_fps();
        let text = format!("FPS: {:.2}", fps);
        draw_text(&text, 10., 50., Self::FONT_SIZE, WHITE);
    }

    fn update(&mut self, dt: f32) {
        self.simulation.update(dt);
    }

    fn draw(&self) {
        clear_background(BLACK);
        self.simulation.draw();
        self.display_stats();
    }
}

#[macroquad::main("Particle simulation")]
async fn main() {
    let mut game_state = GameState::new();

    loop {
        let delta_time = get_frame_time();
        game_state.handle_inputs();
        game_state.update(delta_time);
        game_state.draw();

        next_frame().await;
    }
}
