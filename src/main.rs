use macroquad::prelude::*;

struct GameState {}

impl GameState {
    fn new() -> Self {
        GameState {}
    }

    fn update(&mut self, _dt: f64) {}

    fn draw(&self) {
        clear_background(BLACK);
        draw_circle(0.5, 0.5, 40., WHITE);
    }
}

#[macroquad::main("Particle simulation")]
async fn main() {
    let mut game_state = GameState::new();

    loop {
        let delta_time = get_frame_time() as f64;
        game_state.update(delta_time);
        game_state.draw();
        next_frame().await;
    }
}
