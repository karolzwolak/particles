use macroquad::prelude::*;

#[derive(Clone, Copy)]
/// Positon is not relative to the screen size
/// It is in range [0, 1]
/// (0, 0) is the top left corner
/// We do not need accelaration, is it will be constant and equal to gravity
struct Particle {
    pos: Vec2,
    vel: Vec2,
    color: Color,
}

impl Particle {
    /// The radius of the particle
    /// With a radius of 1, the particle will be as big as the screen
    const RADIUS: f32 = 0.01;

    const GRAVITY: f32 = 0.5;

    fn random_color() -> Color {
        Color::new(
            rand::gen_range(0., 1.),
            rand::gen_range(0., 1.),
            rand::gen_range(0., 1.),
            1.,
        )
    }
    fn new(position: Vec2, velocity: Vec2, color: Option<Color>) -> Self {
        Particle {
            pos: position,
            vel: velocity,
            color: color.unwrap_or_else(|| Self::random_color()),
        }
    }

    fn update(&mut self, dt: f32) {
        self.vel.y += Self::GRAVITY * dt;
        self.pos += self.vel * dt;

        self.clamp_pos_within_borders();
    }

    fn draw(&self) {
        let screen_x = self.pos.x * screen_width();
        let screen_y = self.pos.y * screen_height();

        draw_circle(
            screen_x,
            screen_y,
            Self::RADIUS * screen_width(),
            self.color,
        );
    }

    fn clamp_pos_within_borders(&mut self) {
        self.pos.x = self.pos.x.clamp(0., 1.);
        self.pos.y = self.pos.y.clamp(0., 1.);

        // if particle was colliding with the wall, set vel to zero along this axis
        if self.pos.x == 0. || self.pos.x == 1. {
            self.vel.x = 0.;
        }
        if self.pos.y == 0. || self.pos.y == 1. {
            self.vel.y = 0.;
        }
    }

    fn collides_with(&self, other: &Particle) -> bool {
        let dist = (self.pos - other.pos).length();

        dist < 2. * Self::RADIUS
    }
}

struct Simulation {
    particles: Vec<Particle>,
}

impl Simulation {
    fn new() -> Self {
        Simulation {
            particles: vec![Particle::new(Vec2::new(0., 0.), Vec2::new(0., 0.), None)],
        }
    }

    fn handle_collision(&mut self, a_id: usize, b_id: usize) {
        let a = &self.particles[a_id];
        let b = &self.particles[b_id];

        if a_id == b_id || !a.collides_with(b) {
            return;
        }

        let vec_a_to_b = b.pos - a.pos;

        let a = &mut self.particles[a_id];
        a.vel = Vec2::ZERO;
        a.pos -= vec_a_to_b * 0.5;

        let b = &mut self.particles[b_id];
        b.vel = Vec2::ZERO;
        b.pos += vec_a_to_b * 0.5;
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
            particle.update(dt);
        }
        self.handle_colliosions();
    }

    fn draw(&self) {
        for particle in self.particles.iter() {
            particle.draw();
        }
    }

    fn spawn_particle(&mut self) {
        let pos = Vec2::new(0., 0.1);
        let vel = Vec2::new(0.5, 0.);

        let particle = Particle::new(pos, vel, None);
        self.particles.push(particle);
    }

    fn spawn_at_mouse(&mut self) {
        let mouse = mouse_position();
        let pos = Vec2::new(mouse.0 / screen_width(), mouse.1 / screen_height());
        let vel = Vec2::ZERO;

        let particle = Particle::new(pos, vel, None);
        self.particles.push(particle);
    }
}

struct GameState {
    simulation: Simulation,
}

impl GameState {
    fn new() -> Self {
        GameState {
            simulation: Simulation::new(),
        }
    }

    fn handle_inputs(&mut self) {
        if is_key_pressed(KeyCode::S) {
            self.simulation.spawn_particle();
        }
        if is_mouse_button_pressed(MouseButton::Left) {
            self.simulation.spawn_at_mouse();
        }
    }

    fn update(&mut self, dt: f32) {
        self.simulation.update(dt);
    }

    fn draw(&self) {
        clear_background(BLACK);
        self.simulation.draw();
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
