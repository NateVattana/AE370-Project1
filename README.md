# AE370-Project1

AE 370 group project investigating how passenger count and initial
altitude affect a simplified aircraft model's ability to reach Gimli.

## Requirements

The scripts use Python 3 and Matplotlib.

Install Matplotlib with the same Python interpreter used to run
the scripts:

```bash
python -m pip install matplotlib
```

The code also uses Python's built-in math and os modules.

## Running the code

Open a terminal in the repository folder.

### Main simulation and numerical checks

Run:

```bash
python gimli_simulation.py
```

This script runs:

- The baseline glide simulation.
- The RK4 time-step convergence study.
- The independent projectile-motion validation.
- Initial-altitude, initial-speed, and passenger-count sweeps.
- The 1,000-passenger stress test.
- The fixed-start critical-altitude calculation.

It prints numerical results and saves these figures in the
current working directory:

- baseline_trajectory.png
- rk4_convergence.png
- altitude_sensitivity.png
- airspeed_sensitivity.png
- passenger_sensitivity.png

The script calls plt.show() for the convergence figure. If a
plot window opens and execution pauses, close that window
to allow the remaining studies to run.

### Additional results analysis

Run:

```bash
python results.py
```

This script imports the aircraft model and RK4 solver from
gimli_simulation.py. It calculates baseline arrival conditions,
runs an additional passenger study from 0 to 250 passengers,
compares fixed-start and initially balanced altitude cases,
and evaluates the combined passenger–altitude reachability grid.

It prints numerical results and creates these files inside a
results folder in the current working directory:

- results/A_baseline.png
- results/B_passengers.png
- results/C_altitude.png
- results/D_boundary.png

A_baseline.png is an alternate baseline figure and differs from
baseline_trajectory.png. B_passengers.png uses a different
passenger study from passenger_sensitivity.png. These files
are not interchangeable.

Running either script can overwrite its previously generated
figure files.

## Model scope

The simulations assume zero wind, a prescribed lift coefficient,
and motion along a direct track toward Gimli. They evaluate
whether the modeled aircraft reaches Gimli's horizontal position
before descending to airfield elevation. Approach maneuvers,
touchdown, and braking are not simulated.

The 1,000-passenger case is a deliberately unrealistic
mathematical stress test.
