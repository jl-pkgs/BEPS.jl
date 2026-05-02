#!/usr/bin/env julia
# Usage: julia gen/gen_artifacts.jl [tag]
# Regenerates Artifacts.toml for a new BEPS.c release.
# Requires the release to contain per-platform .tar.gz files:
#   {platform}.tar.gz  e.g. windows-x86_64.tar.gz
# each containing a single shared library file.

using Pkg.Artifacts, Downloads, SHA
using Pkg
using Base.BinaryPlatforms: Platform

const REPO = "CUG-hydro/BEPS.c"
const TAG = get(ARGS, 1, "v2026.05.03")
const BASE_URL = "https://github.com/CUG-hydro/BEPS.c/releases/download/$TAG"

const ARTIFACTS_TOML = joinpath(@__DIR__, "..", "Artifacts.toml")

const PLATFORMS = [
  (platform = "windows-x86_64", arch = "x86_64", os = "windows", lib = "libbeps-windows-x86_64.dll"),
  (platform = "linux-x86_64",   arch = "x86_64", os = "linux",   lib = "libbeps-linux-x86_64.so"),
  (platform = "macos-x86_64",   arch = "x86_64", os = "macos",   lib = "libbeps-macos-x86_64.dylib"),
  (platform = "macos-arm64",    arch = "aarch64", os = "macos",  lib = "libbeps-macos-arm64.dylib"),
]

for p in PLATFORMS
  tarball_url = "$BASE_URL/beps-$(p.platform).tar.gz"
  @info "Downloading" tarball_url

  tmpfile = tempname() * ".tar.gz"
  try
    Downloads.download(tarball_url, tmpfile)
  catch err
    @error "Failed to download" tarball_url err
    continue
  end

  tarball_sha256 = open(io -> bytes2hex(sha256(io)), tmpfile)

  local tree_hash
  try
    tree_hash = create_artifact() do dir
      Pkg.PlatformEngines.unpack(tmpfile, dir)
    end
  catch err
    @error "Failed to unpack" p.platform err
    rm(tmpfile; force=true)
    continue
  end

  rm(tmpfile)

  bind_artifact!(ARTIFACTS_TOML, "libbeps", tree_hash;
    platform = Platform(p.arch, p.os),
    download_info = [(tarball_url, tarball_sha256)],
    lazy = true,
    force = true,
  )

  @info "Bound" p.platform string(tree_hash)
end

@info "Done → Artifacts.toml"
