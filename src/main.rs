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

    grid_last_particle_id: [Option<u32>; Self::GRID_LEN],
    next_particle_neighbor: Vec<Option<u32>>,
}

impl Simulation {
    /// 90% in every direction from center
    const CONSTRAINT_SIZE: f32 = 0.9;

    const CELL_SIZE: f32 = Particle::RADIUS;
    const GRID_ROW_COUNT: usize = (1. / Self::CELL_SIZE) as usize;

    const GRID_LEN: usize = Self::GRID_ROW_COUNT * Self::GRID_ROW_COUNT;

    const SUBSTEPS: u32 = 8;

    const SPAWN_RATE: usize = 10;

    fn dimension() -> f32 {
        screen_width().min(screen_height())
    }

    fn reset_grid(&mut self) {
        for cell in self.grid_last_particle_id.iter_mut() {
            *cell = None;
        }
    }

    fn new() -> Self {
        Simulation {
            particles: Vec::new(),
            grid_last_particle_id: [None; Self::GRID_ROW_COUNT * Self::GRID_ROW_COUNT],
            next_particle_neighbor: Vec::new(),
        }
    }

    fn cell_id(row: usize, col: usize) -> usize {
        row * Self::GRID_ROW_COUNT + col
    }

    fn pos_cell(pos: Vec2) -> (usize, usize) {
        let y = (pos.y + 1.) * 0.5;
        let x = (pos.x + 1.) * 0.5;

        let row = (y * Self::GRID_ROW_COUNT as f32) as usize;
        let col = (x * Self::GRID_ROW_COUNT as f32) as usize;

        (
            row.min(Self::GRID_ROW_COUNT - 1),
            col.min(Self::GRID_ROW_COUNT - 1),
        )
    }

    fn add_particle_to_cell(&mut self, cell_id: usize, particle_id: usize) {
        let last_particle_id = std::mem::replace(
            &mut self.grid_last_particle_id[cell_id],
            Some(particle_id as u32),
        );

        self.next_particle_neighbor[particle_id] = last_particle_id;
    }

    fn populate_grid(&mut self) {
        self.reset_grid();
        for id in 0..self.particles.len() {
            let particle = self.particles[id];
            let (row, col) = Self::pos_cell(particle.pos);
            let cell_id = Self::cell_id(row, col);

            self.add_particle_to_cell(cell_id, id);
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

        let delta_a_b = vec_a_b * (max_dist - dist_a_b) * 0.5 / dist_a_b;

        let a = &mut self.particles[a_id];
        a.pos -= delta_a_b;

        let b = &mut self.particles[b_id];
        b.pos += delta_a_b;
    }

    fn handle_collisions_two_cells(&mut self, cell_id1: usize, cell_id2: usize) {
        let mut next_i = self.grid_last_particle_id[cell_id1];

        while let Some(curr_i) = next_i {
            let mut next_j = self.grid_last_particle_id[cell_id2];
            while let Some(curr_j) = next_j {
                self.handle_collision(curr_i as usize, curr_j as usize);

                next_j = self.next_particle_neighbor[curr_j as usize];
            }
            next_i = self.next_particle_neighbor[curr_i as usize];
        }
    }

    fn handle_collisions(&mut self) {
        let grid_size = Self::GRID_ROW_COUNT as isize;

        for row in 1..Self::GRID_ROW_COUNT - 1 {
            let row_comp = row * Self::GRID_ROW_COUNT;

            for col in 1..Self::GRID_ROW_COUNT - 1 {
                let id1 = row_comp + col;
                if self.grid_last_particle_id[id1].is_none() {
                    continue;
                }
                for row_offset in -1..=1 {
                    let row_comp_offset = row_offset * grid_size;

                    for col_offset in -1..=1 {
                        let offset = row_comp_offset + col_offset;
                        let id2 = (id1 as isize + offset) as usize;

                        if self.grid_last_particle_id[id2].is_none() {
                            continue;
                        }
                        self.handle_collisions_two_cells(id1, id2);
                    }
                }
            }
        }

        for row in 0..Self::GRID_ROW_COUNT {
            let row_comp = row * Self::GRID_ROW_COUNT;

            let mut col = 0;
            while col < Self::GRID_ROW_COUNT {
                if col == 1 && (1..Self::GRID_ROW_COUNT - 1).contains(&row) {
                    col = Self::GRID_ROW_COUNT - 1;
                }
                let id1 = row_comp + col;
                if self.grid_last_particle_id[id1].is_none() {
                    col += 1;
                    continue;
                }
                for row_offset in -1..=1 {
                    let new_row = row as isize + row_offset;
                    if new_row < 0 || new_row >= grid_size {
                        continue;
                    }
                    let row_comp_offset = row_offset * grid_size;

                    for col_offset in -1..=1 {
                        let new_col = col as isize + col_offset;
                        if new_col < 0 || new_col >= grid_size {
                            continue;
                        }
                        let offset = row_comp_offset + col_offset;
                        let id2 = (id1 as isize + offset) as usize;

                        if self.grid_last_particle_id[id2].is_none() {
                            continue;
                        }
                        self.handle_collisions_two_cells(id1, id2);
                    }
                }

                col += 1;
            }
        }
    }

    fn update_once(&mut self, dt: f32) {
        self.populate_grid();
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
        let (row, col) = Self::pos_cell(particle.pos);
        let cell_id = Self::cell_id(row, col);

        self.particles.push(particle);
        self.next_particle_neighbor.push(None);

        self.add_particle_to_cell(cell_id, self.particles.len() - 1);
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
