use rand::Rng;
use std::{thread, time};

//Simulation settings
const NUM_RODS: usize = 5;
const NUM_BOXES: usize = 100;
const GRID_SIZE: usize = 100; //2 * (MOON_RADIUS as usize);

//Graphics settings
const PARTICLE_RADIUS: usize = SCALE / 8;
const WIDTH: usize = 1000; //GRID_SIZE * SCALE;
const HEIGHT: usize = WIDTH; //GRID_SIZE * SCALE;
const SCALE: usize = ((WIDTH as f64) / (GRID_SIZE as f64)) as usize;
const FRAME_TIME: time::Duration = time::Duration::from_millis(16); // ~60 FPS

// Simulation constants - Known
// const MOON_RADIUS: f64 = 142.0;
const SOLAR_WIND_MIN_VELOCITY: f64 = 0.225;

const ROD_DEFLECTION_STRENGTH: f64 = 4.0;
const ROD_MAX_RANGE_SQUARED: f64 = 25.0 * 25.0;
const MIN_ELECTROMAGNETIC_DISTANCE: f64 = 2.0;
const SOLAR_WIND_DEFLECTION_TICK_INTERVAL: u32 = 6;
const SOLAR_WIND_TICK_MULTIPLIER: u32 = 1; //1 when watching, 12 when not?
const COOLDOWN_DISTANCE: f64 = 1.5;
const COOLDOWN_TICKS: f64 = 30.0;

// Changed
const MAX_AGE: f64 = SOLAR_WIND_MIN_VELOCITY * 2.0 * 32.0 * (GRID_SIZE as f64 + 150.0) * 10.0;
const WIND_SPAWN_DISTANCE: f64 = 0.1 * (GRID_SIZE as f64);
const SPAWN_CHANCE: f64 = 1.0; //0.5 * settings.global["cerys-solar-wind-spawn-rate-percentage"].value / 100

#[derive(Debug, Clone, Copy)]
struct Particle {
    x: f64,
    y: f64,
    vx: f64,
    vy: f64,
    age: u32,
    valid: bool,
    irradiation_tick: u32,
    last_checked_inv: u32,
}

#[derive(Debug, Clone, Copy)]
struct Rod {
    gridx: usize,
    gridy: usize,
    x: f64,
    y: f64,
    polarity_fraction: f64,
}

struct Box {
    gridx: usize,
    gridy: usize,
    x: f64,
    y: f64,
    count: u32,
    id: u32,
    last_hit: u32,
}

fn place_rods_and_boxes(rods: &mut Vec<Rod>, boxes: &mut Vec<Box>) {
    let mut rng = rand::rng();

    let mut id = 1;
    while boxes.len() < NUM_BOXES {
        let x = rng.random_range(0..GRID_SIZE - 2);
        let y = rng.random_range(0..GRID_SIZE - 2);

        if !boxes.iter().any(|b| b.gridx == x && b.gridy == y) {
            boxes.push(Box {
                gridx: x,
                gridy: y,
                x: (x as f64) + 0.5,
                y: (y as f64) + 0.5,
                count: 0,
                id: id,
                last_hit: 0,
            });
            id += 1;
        }
    }
    println!("Boxes: {:?}", boxes.len());

    while rods.len() < NUM_RODS {
        let x = rng.random_range(1..GRID_SIZE - 2);
        let y = rng.random_range(1..GRID_SIZE - 2);

        if !rods.iter().any(|r| false) {
            rods.push(Rod {
                gridx: x,
                gridy: y,
                x: (x + 1) as f64,
                y: (y + 1) as f64,
                polarity_fraction: if rng.random_bool(0.5) { 1.0 } else { -1.0 },
            });
        }
    }
}

fn do_tick(tick: u32, particles: &mut Vec<Particle>, rods: &Vec<Rod>, boxes: &mut Vec<Box>) {
    print!("Tick: {} ", tick);
    println!("Particles: {:?}", particles.len());

    let mut rng = rand::rng();

    //Done
    if tick % (1 * SOLAR_WIND_TICK_MULTIPLIER) == 0 {
        tick_1_move_solar_wind(particles)
    }

    if tick % (8 * SOLAR_WIND_TICK_MULTIPLIER) == 0 {
        tick_8_solar_wind_collisions(tick, particles, boxes);
    }

    //Done
    if tick % (SOLAR_WIND_DEFLECTION_TICK_INTERVAL * SOLAR_WIND_TICK_MULTIPLIER) == 0 {
        tick_solar_wind_deflection(particles, rods);
    }

    //Done - 9
    if tick % (9 * SOLAR_WIND_TICK_MULTIPLIER) == 0 {
        if rng.random_range(0.0..1.0) < SPAWN_CHANCE {
            spawn_solar_wind_particle(particles);
        }
    }
    if tick % 240 == 0 {
        tick_240_clean_up_cerys_solar_wind_particles(particles)
    }
}

fn tick_1_move_solar_wind(particles: &mut Vec<Particle>) {
    for p in particles.iter_mut() {
        if !p.valid {
            continue;
        }
        p.x += p.vx;
        p.y += p.vy;
        p.age += 1;
    }
}

fn tick_8_solar_wind_collisions(tick: u32, particles: &mut Vec<Particle>, boxes: &mut Vec<Box>) {
    //this is the big one
    for p in particles.iter_mut() {
        let found_boxes: Vec<&mut Box> = boxes
            .iter_mut()
            .filter(|b| {
                let dx = p.x - b.x;
                let dy = p.y - b.y;
                let d2 = dx * dx + dy * dy;
                d2 < (0.75 * 0.75)
            })
            .collect();
        for b in found_boxes {
            if (!particle_is_in_cooldown(tick, p))
                || ((p.last_checked_inv > 0) && (p.last_checked_inv != b.id))
            {
                p.irradiation_tick = tick;
                p.last_checked_inv = b.id;
                b.count += 1;
                b.last_hit = tick;
                println!(
                    "--------------------------- Box: {:?} Count: {:?}",
                    b.id, b.count
                );
            }
        }
    }
}

fn particle_is_in_cooldown(tick: u32, p: &mut Particle) -> bool {
    if p.irradiation_tick > 0 {
        return false;
    }

    let v2 = p.vx * p.vx + p.vy * p.vy;
    let speed = v2.sqrt();

    let cooldown_time_1 = (COOLDOWN_DISTANCE / speed).round() as u32;
    let cooldown_time_2 = COOLDOWN_TICKS as u32;

    if (tick > p.irradiation_tick + cooldown_time_1)
        || (tick > p.irradiation_tick + cooldown_time_2)
    {
        p.irradiation_tick = 0;
        p.last_checked_inv = 0;
        return false;
    }

    return true;
}

fn tick_solar_wind_deflection(particles: &mut Vec<Particle>, rods: &Vec<Rod>) {
    let mut rng = rand::rng();

    for p in particles.iter_mut() {
        for r in rods.iter() {
            let mut dx = p.x - r.x;
            let mut dy = p.y - r.y;
            let mut d2 = dx * dx + dy * dy;

            if d2 == 0.0 {
                let random_angle = rng.random::<f64>() * 2.0 * std::f64::consts::PI;
                dx = MIN_ELECTROMAGNETIC_DISTANCE * random_angle.cos();
                dy = MIN_ELECTROMAGNETIC_DISTANCE * random_angle.sin();
            } else if d2 < MIN_ELECTROMAGNETIC_DISTANCE * MIN_ELECTROMAGNETIC_DISTANCE {
                let scale = MIN_ELECTROMAGNETIC_DISTANCE / d2.sqrt();
                dx *= scale;
                dy *= scale;
                d2 = dx * dx + dy * dy;
            }

            if (d2 < ROD_MAX_RANGE_SQUARED) && (r.polarity_fraction != 0.0) {
                let deflection = r.polarity_fraction
                    * ROD_DEFLECTION_STRENGTH
                    * (SOLAR_WIND_DEFLECTION_TICK_INTERVAL as f64)
                    / 60.0;
                let dvx = dx / (d2.powf(7.0 / 4.0)) * deflection;
                let dvy = dy / (d2.powf(7.0 / 4.0)) * deflection;
                p.vx += dvx;
                p.vy += dvy;
            }
        }
    }
}

fn spawn_solar_wind_particle(particles: &mut Vec<Particle>) {
    let mut rng = rand::rng();

    for _ in 0..100 {
        let particle = Particle {
            //x: -(WIND_SPAWN_DISTANCE - rng.random_range(0.0..10.0)),
            x: -WIND_SPAWN_DISTANCE * rng.random::<f64>(),
            // y: rng.random_range(0.0..= (2.0 * (MOON_RADIUS + 8.0))),
            y: rng.random_range(0.0..=GRID_SIZE as f64),
            vx: SOLAR_WIND_MIN_VELOCITY + rng.random_range(0.0..1.0) * 0.05,
            vy: 0.3 * (rng.random::<f64>() - 0.5).powi(3),
            age: 0,
            valid: true,
            irradiation_tick: 0,
            last_checked_inv: 0,
        };
        particles.push(particle);
    }
}

fn tick_240_clean_up_cerys_solar_wind_particles(particles: &mut Vec<Particle>) {
    const TOLERANCE: f64 = 1.0;
    for p in particles.iter_mut() {
        if p.age > MAX_AGE as u32 {
            p.valid = false;
        }
        // If they left the sim
        // else if (p.x.abs() > WIND_SPAWN_DISTANCE + 5.0) || (p.y.abs() > WIND_SPAWN_DISTANCE + 5.0)
        else if (p.x < -WIND_SPAWN_DISTANCE)
            || (p.x > (GRID_SIZE as f64) + WIND_SPAWN_DISTANCE)
            || (p.y < -WIND_SPAWN_DISTANCE)
            || (p.y > (GRID_SIZE as f64) + WIND_SPAWN_DISTANCE)
        {
            p.valid = false;
        }
    }

    particles.retain(|p| p.valid);
}

// fn grid_coord_is_on_moon(x: usize, y: usize) -> bool {
//     let dx = (x as f64) - ((GRID_SIZE as f64) / 2.0);
//     let dy = (y as f64) - ((GRID_SIZE as f64) / 2.0);
//     let d2 = dx * dx + dy * dy;
//     d2 < MOON_RADIUS * MOON_RADIUS
// }

fn grid_to_screenx(coord: usize) -> usize {
    coord * SCALE
}