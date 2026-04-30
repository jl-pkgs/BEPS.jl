using BEPS, Test, DataFrames
using BEPS: path_proj, load_forcing, load_lai, load_site_info,
            load_multisite_data, run_multisite
using DataFrames: nrow, names, propertynames

const MULTISITE_DIR = path_proj("examples/multisite")
const SITE_INFO_CSV = joinpath(MULTISITE_DIR, "site_info.csv")

@testset "load_site_info" begin
  si = load_site_info(SITE_INFO_CSV)
  @test si isa DataFrame
  @test nrow(si) == 2
  @test "site_id" in names(si)
  @test "lon"     in names(si)
  @test "weight"  in names(si)        # optional col gets default
  @test si.site_id[1] isa String
  @test si.VegType[1] isa Int
  @test si.lon[1] isa Float64

  # 缺失必要列 → ArgumentError
  si_bad_path = tempname() * ".csv"
  write(si_bad_path, "site_id,lon\nX,1.0\n")
  @test_throws ArgumentError load_site_info(si_bad_path)

  # 文件不存在 → ArgumentError
  @test_throws ArgumentError load_site_info("/no/such/file.csv")
end

@testset "load_forcing" begin
  f = load_forcing(joinpath(MULTISITE_DIR, "CA-Obs", "forcing.txt"))
  @test f isa DataFrame
  @test nrow(f) == 48
  # 列名应已标准化
  @test :Rs   in propertynames(f)
  @test :Tair in propertynames(f)
  @test :rain in propertynames(f)
  @test :Uz   in propertynames(f)
  @test :day  in propertynames(f)
  @test :hour in propertynames(f)
  @test f.Rs[1] isa Float64

  # 文件不存在
  @test_throws ArgumentError load_forcing("/no/such/forcing.txt")
end

@testset "load_lai" begin
  l = load_lai(joinpath(MULTISITE_DIR, "CA-Obs", "lai.txt"))
  @test l isa Vector{Float64}
  @test length(l) == 2
  @test all(l .>= 0)

  # 单行（空格分隔）格式
  tmp = tempname()
  write(tmp, "0.3 0.5 0.7\n")
  l2 = load_lai(tmp)
  @test l2 ≈ [0.3, 0.5, 0.7]

  # 文件不存在
  @test_throws ArgumentError load_lai("/no/such/lai.txt")
end

@testset "load_multisite_data" begin
  si = load_site_info(SITE_INFO_CSV)
  fd, ld = load_multisite_data(si; data_dir=MULTISITE_DIR)

  @test fd isa Dict
  @test ld isa Dict
  @test haskey(fd, "CA-Obs")
  @test haskey(fd, "CA-Qfo")
  @test haskey(ld, "CA-Obs")
  @test fd["CA-Obs"] isa DataFrame
  @test ld["CA-Obs"] isa Vector{Float64}

  # 缺少数据文件的站点被跳过（发出 @warn）
  si_missing = copy(si)
  si_missing.site_id[1] = "NO-SITE"
  fd2, ld2 = load_multisite_data(si_missing; data_dir=MULTISITE_DIR)
  @test !haskey(fd2, "NO-SITE")
end

@testset "run_multisite from path" begin
  # 直接传 site_info.csv 路径
  results = run_multisite(SITE_INFO_CSV; fix_snowpack=false, verbose=false)

  @test results isa Dict
  @test haskey(results, "CA-Obs")
  @test haskey(results, "CA-Qfo")

  for (_, v) in results
    @test v.df_out isa DataFrame
    @test v.df_ET  isa DataFrame
    @test size(v.Tsoil, 2) == 5   # 5 soil layers
    @test size(v.θ,     2) == 5
  end
end

@testset "run_multisite optional columns in site_info" begin
  # 通过 forcing_file / lai_file 指定显式路径
  tmp_csv = tempname() * ".csv"
  obs_path  = joinpath(MULTISITE_DIR, "CA-Obs", "forcing.txt")
  lai_path  = joinpath(MULTISITE_DIR, "CA-Obs", "lai.txt")

  open(tmp_csv, "w") do io
    println(io, "site_id,lon,lat,VegType,SoilType,clumping,Tsoil0,θ0,z_snow0,forcing_file,lai_file")
    println(io, "S1,120.5,30.5,25,8,0.85,2.2,0.41,0.0,$obs_path,$lai_path")
  end

  results = run_multisite(tmp_csv; fix_snowpack=false, verbose=false)
  @test haskey(results, "S1")
end
