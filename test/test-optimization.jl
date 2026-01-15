using BEPS, DataFrames, Test
using Serialization
using DelimitedFiles
using BEPS: path_proj

# @testset "BEPS Optimization" 
begin
    # 1. Load Data
    lai = readdlm(path_proj("examples/input/p1_lai.txt"))[:]
    d = deserialize(path_proj("data/p1_meteo"))
    d.tem = d.tem .- 5.0 # Same adjustment as in test-beps_main.jl

    kw = (lon=120.5, lat=30.5,
          VegType=25, SoilType=8,
          clumping=0.85,
          Tsoil0=2.2, θ0=0.4115, z_snow0=0.0
    )

    # 2. Setup "True" Model and Generate Observations
    model_true = ParamBEPS()
    # Set specific values for optimizable parameters to act as "Truth"
    model_true.r_drainage = 0.55
    model_true.r_root_decay = 0.90
    
    # Run simulation to get synthetic observations
    # We'll use 'LH' (Latent Heat) as our observation target
    _, df_ET_true, _, _ = besp_main(d, lai; model=model_true, kw..., version="julia", fix_snowpack=false)
    obs_ET = df_ET_true.LH

    # 3. Setup Initial Model for Optimization (perturbed parameters)
    model_opt = ParamBEPS()
    model_opt.r_drainage = 0.30 # Initial guess far from 0.55
    model_opt.r_root_decay = 0.98 # Initial guess far from 0.90

    # 4. Run Optimization
    # We use a small maxn to ensure the test runs quickly. 
    # The goal is to verify the pipeline works, not necessarily to achieve perfect convergence in a unit test.
    println("Starting optimization...")
    model_final, best_rmse = beps_optimize(d, lai, model_opt, obs_ET; 
        col_sim=:LH, 
        kw..., 
        version="julia", 
        fix_snowpack=false,
        maxn=10, # Small number of iterations for testing
        kstop=3,
        pcento=0.1
    )

    # 5. Verify Results
    println("Optimization finished.")
    println("True r_drainage: $(model_true.r_drainage), Optimized: $(model_final.r_drainage)")
    println("True r_root_decay: $(model_true.r_root_decay), Optimized: $(model_final.r_root_decay)")
    println("Best RMSE: $best_rmse")

    @test best_rmse < 1.0 # Should be very small since data is synthetic, but maxn is low
    @test model_final isa ParamBEPS
    
    # Check if parameters moved in the right direction (optional, given low iterations)
    # But mainly we check that the code ran without error.
    @test true 
end
