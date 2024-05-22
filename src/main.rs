use macroquad::prelude::*;

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
    const RADIUS: f32 = 0.002;

    const GRAVITY: f32 = 0.5;
    const GRAVITY_VEC: Vec2 = Vec2::new(0., Self::GRAVITY);
    const SCALE: f32 = 0.5 * Simulation::CONSTRAINT_SIZE;

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
            prev_pos: position - initial_velocity * get_frame_time() / Simulation::SUBSTEPS as f32,
            color: color.unwrap_or_else(Self::random_color),
        }
    }

    fn update_verlet(&mut self, dt: f32) {
        let new_pos = 2. * self.pos - self.prev_pos + Self::GRAVITY_VEC * dt * dt;
        self.prev_pos = self.pos;
        self.pos = new_pos;
    }

    fn draw(&self) {
        // we use half of the dimension, because the position is in range [-1, 1]
        let half_dim = Simulation::dimension() * Self::SCALE;

        let screen_x = screen_width() * 0.5 + self.pos.x * half_dim;
        let screen_y = screen_height() * 0.5 + self.pos.y * half_dim;

        draw_circle(screen_x, screen_y, Self::RADIUS * half_dim, self.color);
    }

    fn constraint(&mut self) {
        self.pos.x = self.pos.x.clamp(-1., 1.);
        self.pos.y = self.pos.y.clamp(-1., 1.);
    }
}

struct Simulation {
    particles: Vec<Particle>,
}

impl Simulation {
    /// 90% in every direction from center
    const CONSTRAINT_SIZE: f32 = 0.9;

    const SUBSTEPS: u32 = 8;

    const SPAWN_RATE: usize = 10;

    fn dimension() -> f32 {
        screen_width().min(screen_height())
    }

    fn new() -> Self {
        Simulation {
            particles: Vec::new(),
        }
    }

    fn handle_collision(&mut self, a_id: usize, b_id: usize) {
        if a_id == b_id {
            return;
        }
        let a = &self.particles[a_id];
        let b = &self.particles[b_id];

        let vec_a_b = b.pos - a.pos;
        let dist_a_b = vec_a_b.length();

        let max_dist = 2. * Particle::RADIUS;

        if dist_a_b >= max_dist {
            return;
        }
        let min_dist = 0.01 * Particle::RADIUS;
        let dist_a_b = dist_a_b.max(min_dist);

        let delta_a_b = vec_a_b * (max_dist - dist_a_b) * 0.5 / dist_a_b;

        let a = &mut self.particles[a_id];
        a.pos -= delta_a_b;

        let b = &mut self.particles[b_id];
        b.pos += delta_a_b;
    }

    fn hande_collisions_between_groups(
        &mut self,
        split: Vec2,
        split_along_x: bool,
        a: &[(usize, Vec2)],
        b: &[(usize, Vec2)],
    ) {
        let a_close = a
            .iter()
            .rev()
            .take_while(|(_, pos)| {
                let diff = *pos - split;
                let val = if split_along_x { diff.x } else { diff.y };
                val <= Particle::RADIUS
            })
            .map(|(id, _)| *id)
            .collect::<Vec<_>>();

        let b_close = b
            .iter()
            .take_while(|(_, pos)| {
                let diff = split - *pos;
                let val = if split_along_x { diff.x } else { diff.y };
                val <= Particle::RADIUS
            })
            .map(|(id, _)| *id)
            .collect::<Vec<_>>();

        for a_id in a_close.iter() {
            for b_id in b_close.iter() {
                self.handle_collision(*a_id, *b_id);
            }
        }
    }

    fn divide_particles<'a>(
        &self,
        split_along_x: bool,
        particles: &'a mut [(usize, Vec2)],
    ) -> (Vec2, &'a mut[(usize, Vec2)], &'a mut [(usize, Vec2,)]) {
        particles.sort_unstable_by(|a, b| {
            let pos_a = a.1;
            let pos_b = b.1;

            if split_along_x {
                pos_a.y.partial_cmp(&pos_b.y).unwrap()
            } else {
                pos_a.x.partial_cmp(&pos_b.x).unwrap()
            }
        });

        let split_at = particles.len() / 2;
        let split_pos = if split_at % 2 == 0 {
            particles[split_at].1
        } else {
            (particles[split_at - 1].1 + particles[split_at].1) * 0.5
        };

        let (a, b) = particles.split_at_mut(split_at);

        (split_pos, a, b)
    }

    fn divide_handle_collision(&mut self, min: Vec2, max: Vec2, particles: &mut [(usize, Vec2)]) {
        if particles.len() <= 1 {
            return;
        }
        let diff = max - min;
        let split_along_x = diff.x > diff.y;

        let (split, a, b) = self.divide_particles(split_along_x, particles);
        let mut split_min = min;
        let mut split_max = max;

        if split_along_x {
            split_min.y = split.y;
            split_max.y = split.y;
        } else {
            split_min.x = split.x;
            split_max.x = split.x;
        }

        self.hande_collisions_between_groups(split, split_along_x, &a, &b);

        self.divide_handle_collision(min, split_max, a);
        self.divide_handle_collision(split_min, max, b);
    }

    fn handle_collisions(&mut self) {
        let min = Vec2::new(-1., -1.);
        let max = Vec2::new(1., 1.);
        let mut particles = self
            .particles
            .iter()
            .enumerate()
            .map(|(i, p)| (i, p.pos))
            .collect::<Vec<_>>();

        self.divide_handle_collision(min, max, &mut particles);
    }

    fn update_once(&mut self, dt: f32) {
        self.handle_collisions();
        for particle in self.particles.iter_mut() {
            particle.update_verlet(dt);
            particle.constraint();
        }
    }

    fn update_substeps(&mut self, dt: f32) {
        let substep_dt = dt / Self::SUBSTEPS as f32;
        for _ in 0..Self::SUBSTEPS {
            self.update_once(substep_dt);
        }
    }

    fn draw(&self) {
        for particle in self.particles.iter() {
            particle.draw();
        }
    }

    fn add_particle(&mut self, particle: Particle) {
        self.particles.push(particle);
    }

    fn spawn_particles(&mut self) {
        let vel = Vec2::new(1.0, 0.);
        let pos = Vec2::new(-0.9, -0.6);

        for i in 0..Self::SPAWN_RATE {
            let y_offset = i as f32 * Particle::RADIUS * 2.0;
            let pos = pos + Vec2::new(0., y_offset);

            let particle = Particle::new(pos, vel, None);
            self.add_particle(particle);
        }
    }

    fn spawn_at_mouse(&mut self) {
        let mouse = mouse_position();
        let dimension = Self::dimension();

        let x = (mouse.0 - screen_width() * 0.5) / dimension;
        let y = (mouse.1 - screen_height() * 0.5) / dimension;
        let pos = Vec2::new(x, y);
        let vel = Vec2::ZERO;

        let particle = Particle::new(pos, vel, None);
        self.add_particle(particle);
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
            self.simulation.spawn_particles();
        }
        if is_key_down(KeyCode::A) {
            self.simulation.spawn_particles();
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
        self.simulation.update_substeps(dt);
    }

    fn draw(&self) {
        clear_background(BLACK);
        self.simulation.draw();
        self.display_stats();
    }
}

fn config() -> Conf {
    Conf {
        window_title: "Particle simulation".to_string(),
        fullscreen: true,
        ..Default::default()
    }
}

#[macroquad::main(config)]
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
