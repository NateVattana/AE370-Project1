import math

#aircraft parameters
m_fixed = 80000        #fixed aircraft mass, including crew (kg)
N_p = 61               #number of passengers
m_p = 100              #mass per passenger, including baggage (kg)
m = m_fixed + N_p*m_p  #total aircraft mass (kg)

S = 283.3              #wing area (m^2)
C_L = 0.26             #constant lift coefficient
C_D0 = 0.020           #zero-lift drag coefficient
K = 0.0498             #induced-drag factor

#atmosphere parameters
rho_SL = 1.225         #sea-level reference density (kg/m^3)
H = 8500               #atmospheric scale height (m)
g = 9.81               #gravitational acceleration (m/s^2)

#initial conditions
x0 = 0                 #init horizontal position (m)
h0 = 10668             #init altitude (m)
V0 = 246.77            #init airspeed (m/s)
gamma0 = 0             #init flight-path angle (radians)

#airport and wind conditions
x_G = 72420            #horizontal distance to Gimli (m)
h_G = 230              #gimli airfield elevation (m)
w = 0                  #horizontal wind velocity (m/s)


#governing diff eqs
def aircraft_dynamics(t, y, mass=m):
    #state variables
    x, h, V, gamma = y

    #atmospheric density (varies with altitude)
    rho = rho_SL * math.exp(-h / H)

    #aerodynamic lift/drag
    L = 0.5 * rho * V**2 * S * C_L
    D = 0.5 * rho * V**2 * S * (C_D0 + K * C_L**2)

    #four governing diff eqs
    dxdt = V * math.cos(gamma) + w                         #horizontal v
    dhdt = V * math.sin(gamma)                             #vertical v
    dVdt = -D / mass - g * math.sin(gamma)                 #ROC of airspeed
    dgammadt = L / (mass * V) - g * math.cos(gamma) / V    #ROC of flight path angle

    return [dxdt, dhdt, dVdt, dgammadt]


#fourth-order Runge-Kutta (RK4)
def rk4_step(t, y, dt, f=aircraft_dynamics):
    #first slope: beginning of the time step
    k1 = f(t, y)

    #second slope: first midpoint estimate
    y2 = [y[i] + (dt / 2) * k1[i] for i in range(4)]
    k2 = f(t + dt / 2, y2)

    #third slope: second midpoint estimate
    y3 = [y[i] + (dt / 2) * k2[i] for i in range(4)]
    k3 = f(t + dt / 2, y3)

    #fourth slope: end of the time step
    y4 = [y[i] + dt * k3[i] for i in range(4)]
    k4 = f(t + dt, y4)

    #combine slopes to calculate the next state
    y_next = [
        y[i] + (dt / 6) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])
        for i in range(4)
    ]

    return y_next


#full glide simulation w adjustable paramters
def simulate_glide(dt=1.0, initial_altitude=h0,
                   initial_speed=V0, passengers=N_p):

    #calculate aircraft mass for this scenario
    scenario_mass = m_fixed + passengers * m_p

    #use the selected mass in the governing equations
    def scenario_dynamics(t, y):
        return aircraft_dynamics(t, y, mass=scenario_mass)

    #initial conds
    t = 0.0
    y = [x0, initial_altitude, initial_speed, gamma0]

    #store flight history
    times = [t]
    states = [y.copy()]

    max_time = 7200.0

    #advance until reaching airfield elevation
    while y[1] > h_G and t < max_time:

        if y[2] <= 0:
            raise ValueError("Airspeed became zero or negative.")

        step = min(dt, max_time - t)

        #RK4 using the scenario-specific aircraft mass
        y_next = rk4_step(t, y, step, f=scenario_dynamics)

        #locate airfield elevation crossing
        if y_next[1] <= h_G:
            low = 0.0
            high = step

            for _ in range(40):
                mid = (low + high) / 2

                y_mid = rk4_step(t, y, mid, f=scenario_dynamics)

                if y_mid[1] > h_G:
                    low = mid
                else:
                    high = mid

            step = high
            y_next = rk4_step(t, y, step, f=scenario_dynamics)

            t += step
            y = y_next
            times.append(t)
            states.append(y.copy())
            break

        t += step
        y = y_next

        times.append(t)
        states.append(y.copy())

    else:
        raise RuntimeError(
            "Simulation ended before reaching airfield elevation."
        )

    return times, states


#run simulation
times, states = simulate_glide(dt=1.0)

#final conditions
final_state = states[-1]
final_x = final_state[0]
final_h = final_state[1]
final_V = final_state[2]
final_gamma = final_state[3]

#compare glide distance with distance to Gimli
distance_margin = final_x - x_G

print("\n--- Glide Results ---")
print(f"Flight time: {times[-1]:.2f} seconds")
print(f"Horizontal distance: {final_x:.2f} m")
print(f"Final altitude: {final_h:.2f} m")
print(f"Final airspeed: {final_V:.2f} m/s")
print(f"Final flight-path angle: {math.degrees(final_gamma):.2f} degrees")
print(f"Distance margin: {distance_margin:.2f} m")

if distance_margin >= 0:
    print("Result: Gimli was within the modeled glide range.")
else:
    print("Result: Aircraft fell short of Gimli.")


#convergence
#time steps to investigate
dt_values = [2.0, 1.0, 0.5, 0.25, 0.125]

#compute reference solution
dt_reference = 0.01
_, reference_states = simulate_glide(dt=dt_reference)
x_reference = reference_states[-1][0]

print("\n--- RK4 Convergence ---")
print(f"Reference distance: {x_reference:.8f} m")
print(f"{'dt (s)':<12}{'Distance (m)':<20}{'Error (m)':<20}{'Ratio':<12}")

previous_error = None

for dt_test in dt_values:

    _, test_states = simulate_glide(dt=dt_test)

    x_test = test_states[-1][0]

    #error relative to reference solution
    error = abs(x_test - x_reference)

    #error reduction compared with previous time step
    if previous_error is None:
        ratio = None
    elif error > 0:
        ratio = previous_error / error
    else:
        ratio = None

    ratio_text = f"{ratio:.3f}" if ratio is not None else "---"

    print(f"{dt_test:<12.3f}{x_test:<20.8f}{error:<20.10f}{ratio_text:<12}")

    previous_error = error


#independent validation: ideal projectile motion
#remove aerodynamic lift and drag for this test only

def projectile_dynamics(t, y):
    x, h, V, gamma = y

    dxdt = V * math.cos(gamma)
    dhdt = V * math.sin(gamma)
    dVdt = -g * math.sin(gamma)
    dgammadt = -g * math.cos(gamma) / V

    return [dxdt, dhdt, dVdt, dgammadt]


#validation initial conditions
t = 0.0
dt_test = 0.1
t_final = 100.0
y = [x0, h0, V0, gamma0]

#integrate simplified system
while t < t_final - 1e-12:
    step = min(dt_test, t_final - t)
    y = rk4_step(t, y, step, f=projectile_dynamics)
    t += step

#exact analytical solution at t = 100 s
x_exact = x0 + V0 * math.cos(gamma0) * t_final
h_exact = h0 + V0 * math.sin(gamma0) * t_final - 0.5 * g * t_final**2

#errors
error_x = abs(y[0] - x_exact)
error_h = abs(y[1] - h_exact)

print("\n--- Independent Validation ---")
print(f"Exact horizontal position: {x_exact:.6f} m")
print(f"RK4 horizontal position:   {y[0]:.6f} m")
print(f"Horizontal error:          {error_x:.8f} m")

print(f"Exact altitude:            {h_exact:.6f} m")
print(f"RK4 altitude:              {y[1]:.6f} m")
print(f"Altitude error:            {error_h:.8f} m")



#RK4 convergence plot
import matplotlib.pyplot as plt

plot_dt = []
plot_errors = []

#numerical errors for each time step
for dt_test in dt_values:
    _, test_states = simulate_glide(dt=dt_test)
    x_test = test_states[-1][0]

    plot_dt.append(dt_test)
    plot_errors.append(abs(x_test - x_reference))

#plot numerical error
plt.figure(figsize=(7, 5))

plt.loglog(
    plot_dt,
    plot_errors,
    "o-",
    label="RK4 numerical error"
)

#theoretical fourth-order convergence
reference_curve = [
    plot_errors[0] * (step / plot_dt[0])**4
    for step in plot_dt
]

plt.loglog(
    plot_dt,
    reference_curve,
    "--",
    label="Fourth-order reference"
)

plt.xlabel("Time step (s)")
plt.ylabel("Absolute glide-range error (m)")
plt.title("RK4 Convergence Study - Gimli Glider")
plt.legend()
plt.grid(True, which="both", alpha=0.3)

plt.tight_layout()
plt.savefig("rk4_convergence.png", dpi=300)
plt.show()


#test adjustable simulation parameters
_, baseline = simulate_glide(dt=1.0)

_, lower_altitude = simulate_glide(
    dt=1.0,
    initial_altitude=h0 * 0.95
)

_, higher_mass = simulate_glide(
    dt=1.0,
    passengers=92
)

print("\n--- Adjustable Simulation Test ---")
print(f"Baseline range: {baseline[-1][0]:.2f} m")
print(f"5% lower altitude range: {lower_altitude[-1][0]:.2f} m")
print(f"92 passenger range: {higher_mass[-1][0]:.2f} m")



#automated parameter sweeps
#test how initial altitude, airspeed, and passenger count affect the aircraft's maximum horizontal glide range.

#percentages of the baseline values
percentages = [0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20]

#store results for later analysis
altitude_results = []
speed_results = []
passenger_results = []

#1. initial altitude sweep
print("\n--- Initial Altitude Sweep ---")

for factor in percentages:
    test_altitude = h0 * factor

    _, states_test = simulate_glide(
        dt=1.0,
        initial_altitude=test_altitude
    )

    glide_range = states_test[-1][0]
    altitude_results.append((factor, glide_range))

    print(f"{factor*100:.0f}% altitude: {glide_range:.2f} m")


#2. initial airspeed sweep
print("\n--- Initial Airspeed Sweep ---")

for factor in percentages:
    test_speed = V0 * factor

    _, states_test = simulate_glide(
        dt=1.0,
        initial_speed=test_speed
    )

    glide_range = states_test[-1][0]
    speed_results.append((factor, glide_range))

    print(f"{factor*100:.0f}% speed: {glide_range:.2f} m")


#3. passenger count sweep
print("\n--- Passenger Count Sweep ---")

for factor in percentages:
    test_passengers = round(N_p * factor)

    _, states_test = simulate_glide(
        dt=1.0,
        passengers=test_passengers
    )

    glide_range = states_test[-1][0]
    passenger_results.append((factor, glide_range))

    print(
        f"{test_passengers} passengers: "
        f"{glide_range:.2f} m"
    )



# Find the critical initial altitude using bisection

def find_critical_altitude(dt=1.0, tolerance=0.01):

    # Lower and upper altitude bounds (meters)
    low = 3000.0
    high = h0

    # Check that the bounds have opposite outcomes
    _, low_states = simulate_glide(
        dt=dt, initial_altitude=low
    )
    _, high_states = simulate_glide(
        dt=dt, initial_altitude=high
    )

    low_margin = low_states[-1][0] - x_G
    high_margin = high_states[-1][0] - x_G

    if low_margin >= 0 or high_margin < 0:
        raise ValueError(
            "Altitude bounds do not bracket the failure threshold."
        )

    # Narrow down the critical altitude
    for _ in range(50):
        mid = (low + high) / 2

        _, mid_states = simulate_glide(
            dt=dt, initial_altitude=mid
        )

        margin = mid_states[-1][0] - x_G

        if abs(margin) < tolerance:
            return mid, margin

        if margin >= 0:
            high = mid
        else:
            low = mid

    # Return the best estimate
    critical_altitude = (low + high) / 2

    _, final_states = simulate_glide(
        dt=dt, initial_altitude=critical_altitude
    )

    final_margin = final_states[-1][0] - x_G

    return critical_altitude, final_margin


#calculate altitude margin for error
critical_altitude, critical_margin = find_critical_altitude()

altitude_reduction = (
    (h0 - critical_altitude) / h0
) * 100

print("\n--- Critical Initial Altitude ---")
print(f"Critical altitude: {critical_altitude:.2f} m")
print(f"Critical altitude: {critical_altitude / 0.3048:.2f} ft")
print(f"Allowable altitude reduction: {altitude_reduction:.2f}%")
print(f"Distance error at threshold: {critical_margin:.4f} m")


#parameter sensitivity graphs
#uses the parameter sweep results already calculated

def plot_sensitivity(results, xlabel, title, filename,
                     critical_value=None):

    #separate parameter values and glide distances
    x_values = [item[0] * 100 for item in results]
    ranges_km = [item[1] / 1000 for item in results]

    plt.figure(figsize=(7, 5))

    #predicted glide range
    plt.plot(
        x_values, ranges_km, "o-",
        label="Predicted glide range"
    )

    #distance required to reach Gimli
    plt.axhline(
        y=x_G / 1000,
        color="red",
        linestyle="--",
        label="Distance to Gimli"
    )

    #critical altitude marker, if applicable
    if critical_value is not None:
        plt.scatter(
            [critical_value],
            [x_G / 1000],
            color="black",
            zorder=5,
            label="Critical threshold"
        )

    plt.xlabel(xlabel)
    plt.ylabel("Maximum modeled glide range (km)")
    plt.title(title)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    plt.savefig(filename, dpi=300)
    plt.close()


#1. initial altitude sensitivity
critical_percent = critical_altitude / h0 * 100

plot_sensitivity(
    altitude_results,
    "Initial altitude (% of baseline)",
    "Effect of Initial Altitude on Glide Range",
    "altitude_sensitivity.png",
    critical_value=critical_percent
)

#2. initial airspeed sensitivity
plot_sensitivity(
    speed_results,
    "Initial airspeed (% of baseline)",
    "Effect of Initial Airspeed on Glide Range",
    "airspeed_sensitivity.png"
)

#3. passenger count sensitivity
plot_sensitivity(
    passenger_results,
    "Passenger count (% of baseline)",
    "Effect of Passenger Count on Glide Range",
    "passenger_sensitivity.png"
)

print("\nSensitivity graphs saved successfully.")
