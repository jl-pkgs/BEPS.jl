cd src
npx repomix . --style markdown --include '*.jl,*/*.jl,*/*/*.jl' --ignore 'backup, clang' -o ../BEPS_v0.1.9.md
cd ..
