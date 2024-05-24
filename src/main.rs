use std::collections::HashSet;

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
        let half_dim = Simulation::half_dim();
        let screen_pos = Simulation::to_screen_pos(self.pos);

        draw_circle(
            screen_pos.x,
            screen_pos.y,
            Self::RADIUS * half_dim,
            self.color,
        );
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
    const SCALE: f32 = 0.5 * Self::CONSTRAINT_SIZE;

    const SUBSTEPS: u32 = 8;

    const SPAWN_RATE: usize = 10;

    fn dimension() -> f32 {
        screen_width().min(screen_height())
    }

    fn half_dim() -> f32 {
        Simulation::dimension() * Self::SCALE
    }

    fn to_screen_pos(local_pos: Vec2) -> Vec2 {
        let half_dim = Self::half_dim();

        let screen_x = screen_width() * 0.5 + local_pos.x * half_dim;
        let screen_y = screen_height() * 0.5 + local_pos.y * half_dim;

        Vec2::new(screen_x, screen_y)
    }

    fn to_local_pos(screen: Vec2) -> Vec2 {
        let half_dim = Self::half_dim();

        let x = screen.x - screen_width() * 0.5;
        let y = screen.y - screen_height() * 0.5;

        Vec2::new(x, y) / half_dim
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
        let vec_a_b = self.particles[b_id].pos - self.particles[a_id].pos;
        let dist_a_b = vec_a_b.length();

        let max_dist = 2. * Particle::RADIUS;

        if dist_a_b >= max_dist {
            return;
        }
        let min_dist = 0.01 * Particle::RADIUS;
        let dist_a_b = dist_a_b.max(min_dist);

        let delta_a_b = vec_a_b * (max_dist - dist_a_b) * 0.5 / dist_a_b;

        self.particles[a_id].pos -= delta_a_b;
        self.particles[b_id].pos += delta_a_b;
    }

    fn hande_collisions_between_groups(
        &mut self,
        split: Vec2,
        split_along_x: bool,
        a: &[usize],
        b: &[usize],
    ) {
        let within_range = |p_id: &usize| {
            let diff = self.particles[*p_id].pos - split;
            let val = if split_along_x { diff.y } else { diff.x };
            val.abs() <= 2. * Particle::RADIUS
        };
        let a_close = a
            .iter()
            .copied()
            .rev()
            .take_while(within_range)
            .collect::<Vec<_>>();

        let b_close = b
            .iter()
            .copied()
            .take_while(within_range)
            .collect::<Vec<_>>();

        if a_close.is_empty() || b_close.is_empty() {
            return;
        }
        for a in a_close.iter() {
            for b in b_close.iter() {
                self.handle_collision(*a, *b);
            }
        }
    }

    fn divide_particles<'a>(
        &self,
        split_along_x: bool,
        particles_x_sorted: &'a mut [usize],
        particles_y_sorted: &'a mut [usize],
    ) -> (
        Vec2,
        (&'a mut [usize], &'a mut [usize]),
        (&'a mut [usize], &'a mut [usize]),
    ) {
        let split_at = particles_x_sorted.len() / 2;

        let (primary_particles, secondary_particles) = if split_along_x {
            (particles_y_sorted, particles_x_sorted)
        } else {
            (particles_x_sorted, particles_y_sorted)
        };

        let split_pos = if primary_particles.len() % 2 == 1 {
            let id = primary_particles[split_at];
            self.particles[id].pos
        } else {
            let id = primary_particles[split_at];
            let prev_id = primary_particles[split_at - 1];
            (self.particles[id].pos + self.particles[prev_id].pos) * 0.5
        };

        let primary_div = primary_particles.split_at_mut(split_at);
        let primary_a_indeces = primary_div.0.iter().copied().collect::<HashSet<_>>();

        let secondary_clone = secondary_particles.iter().copied().collect::<Vec<_>>();
        let secondary_div = secondary_particles.split_at_mut(split_at);

        let mut a_count = 0;
        let mut b_count = 0;
        for p_id in secondary_clone {
            // additional condition to skip cheking the set if a is full
            if a_count < split_at && primary_a_indeces.contains(&p_id) {
                secondary_div.0[a_count] = p_id;
                a_count += 1;
            } else {
                secondary_div.1[b_count] = p_id;
                b_count += 1;
            }
        }

        let (x_sorted_div, y_sorted_div) = if split_along_x {
            (secondary_div, primary_div)
        } else {
            (primary_div, secondary_div)
        };

        (split_pos, x_sorted_div, y_sorted_div)
    }

    fn divide_handle_collision(
        &mut self,
        min: Vec2,
        max: Vec2,
        particles_x_sorted: &mut [usize],
        particles_y_sorted: &mut [usize],
    ) {
        if particles_x_sorted.len() <= 1 {
            return;
        }

        let diff = max - min;
        let split_along_x = diff.x < diff.y;

        let (split, x_sorted_div, y_sorted_div) =
            self.divide_particles(split_along_x, particles_x_sorted, particles_y_sorted);

        let mut split_min = min;
        let mut split_max = max;

        if split_along_x {
            split_min.y = split.y;
            split_max.y = split.y;
        } else {
            split_min.x = split.x;
            split_max.x = split.x;
        }

        let primary_div = if split_along_x {
            &y_sorted_div
        } else {
            &x_sorted_div
        };

        self.hande_collisions_between_groups(split, split_along_x, primary_div.0, primary_div.1);

        self.divide_handle_collision(min, split_max, x_sorted_div.0, y_sorted_div.0);
        self.divide_handle_collision(split_min, max, x_sorted_div.1, y_sorted_div.1);
    }

    fn handle_collisions(&mut self) {
        let min = Vec2::new(-1., -1.);
        let max = Vec2::new(1., 1.);

        let indeces = (0..self.particles.len()).collect::<Vec<_>>();

        let mut particles_x_sorted = indeces.clone();
        let mut particles_y_sorted = indeces;

        particles_x_sorted.sort_unstable_by(|a_id, b_id| {
            let pos_a = self.particles[*a_id].pos;
            let pos_b = self.particles[*b_id].pos;

            pos_a.x.partial_cmp(&pos_b.x).unwrap()
        });
        particles_y_sorted.sort_unstable_by(|a_id, b_id| {
            let pos_a = self.particles[*a_id].pos;
            let pos_b = self.particles[*b_id].pos;

            pos_a.y.partial_cmp(&pos_b.y).unwrap()
        });

        self.divide_handle_collision(min, max, &mut particles_x_sorted, &mut particles_y_sorted);
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
        let mouse = Vec2::from(mouse_position());

        let pos = Self::to_local_pos(mouse);
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
