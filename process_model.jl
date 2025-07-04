
using ModelingToolkit
using DifferentialEquations

struct ModelSymbols
    parameters::Vector{}
    variables::Vector{}
    controls::Vector{}
end


struct ProcessModel
    n_steps::Int
    symbols::ModelSymbols
    model::ODESystem
end


function create_process_model(n_steps::Int)
    @independent_variables t
    @variables y1(t) y2(t) y3(t) y4(t) y5(t) y6(t) y7(t) y8(t) y9(t)
    @variables G(t) B(t) C(t) Cs(t) Cm(t) T(t)
    @variables C_high(t) C_low(t) C_high0(t) dT_pos(t)  # over and under concentration variables

    # Control discretization parameters
    @parameters ρ [tunable = false] kv [tunable = false] kb [tunable = false] kg [tunable = false] Eg_R [tunable = false] Eb_R [tunable = false] b [tunable = false] g [tunable = false] h_offset [tunable = false]

    # Generate control parameter symbols dynamically
    u_param_symbols = [Symbol("u_$i") for i in 1:n_steps+1]

    # Create parameter variables
    eval(:(@parameters $(u_param_symbols...)))

    # Bound parameters
    D = Differential(t)

    # Helper function to create piecewise control function using recursion
    function create_piecewise_control(control_params, dt_val, t_var)
        n_pieces = length(control_params)

        function build_recursive(i)
            if i == n_pieces
                return control_params[i]
            else
                return ifelse(t_var < i * dt_val, control_params[i] + (control_params[i+1] - control_params[i]) * (t - (i - 1) * dt_val) / dt_val, build_recursive(i + 1))
                #return ifelse(t_var < i * dt_val, control_params[i], build_recursive(i + 1))
            end
        end

        return build_recursive(1)
    end

    function create_check_pos(control_params, dt_val, t_var)
        n_pieces = length(control_params)

        function build_recursive(i)
            if i == n_pieces
                return 0.0
            else
                return ifelse(t_var < i * dt_val, ifelse(control_params[i+1] > control_params[i], control_params[i+1] - control_params[i], 0), build_recursive(i + 1))
            end
        end

        return build_recursive(1)
    end


    # Control variable definitions using recursive function
    u_params = [eval(sym) for sym in u_param_symbols]
    dt = 1800 / n_steps
    T_eq = T ~ create_piecewise_control(u_params, dt, t)

    # Define saturation and metastable concentrations as functions of temperature
    Cs_eq = Cs ~ 6.29e-2 + 2.46e-3 * (T - 273.15) - 7.14e-6 * (T - 273.15)^2
    Cm_eq = Cm ~ 7.76e-2 + 2.46e-3 * (T - 273.15) - 8.10e-6 * (T - 273.15)^2

    # Concentration (C is y1)
    C_eq = C ~ y1
    # Supersaturation ratio based on saturation concentration
    supersaturation_ratio = ifelse(C > Cs, (C - Cs) / Cs, 0.0)

    B_eq = B ~ kb * exp(-Eb_R / (T)) * supersaturation_ratio^b * (y5 + y9)

    G_eq = G ~ kg * exp(-Eg_R / (T)) * supersaturation_ratio^g


    eqs = [
        # Solute Concentration
        D(y1) ~ -3 * ρ * kv * G * (y4 + y8),

        # Zeroth Moment of Seeded Crystals
        D(y2) ~ 0,

        # First Moment of Seeded Crystals
        D(y3) ~ G * y2,

        # Second Moment of Seeded Crystals
        D(y4) ~ 2 * G * y3,

        # Third Moment of Seeded Crystals
        D(y5) ~ 3 * G * y4,

        # Zeroth Moment of Nucleated Crystals
        D(y6) ~ B,

        # First Moment of Nucleated Crystals
        D(y7) ~ G * y6,

        # Second Moment of Nucleated Crystals
        D(y8) ~ 2 * G * y7,

        # Third Moment of Nucleated Crystals
        D(y9) ~ 3 * G * y8,

        # bound violation tracking equations
        D(C_low) ~ ifelse((Cs - C) > 0, (Cs - C), 0),
        D(C_high) ~ ifelse((C + h_offset - Cm) > 0, (C + h_offset - Cm), 0),
        D(C_high0) ~ ifelse((C - Cm) > 0, (C - Cm), 0),
        D(dT_pos) ~ create_check_pos(u_params, dt, t) * 100,  # Check if temperature is increasing
        Cs_eq,
        Cm_eq,
        C_eq,
        B_eq,
        G_eq,
        T_eq
    ]

    variables_list = [y1, y2, y3, y4, y5, y6, y7, y8, y9, G, B, C, Cs, Cm, T, C_high, C_low, C_high0, dT_pos]

    parameters_list = [ρ, kv, kb, kg, Eg_R, Eb_R, b, g, h_offset]

    all_params = [parameters_list..., u_params...]

    # Create the ODESystem
    @named crystallization_model = ODESystem(eqs, t, variables_list, all_params)

    # Simplify the system
    crystallization_model = structural_simplify(crystallization_model)

    # Create control_symbols dictionary
    control_list = u_params

    return crystallization_model, parameters_list, variables_list, control_list
end

function ProcessModel(n_steps::Int)
    crystallization_model, parameters_list, variables_list, control_list = create_process_model(n_steps)
    symbols = ModelSymbols(parameters_list, variables_list, control_list)
    return ProcessModel(n_steps, symbols, crystallization_model)
end
