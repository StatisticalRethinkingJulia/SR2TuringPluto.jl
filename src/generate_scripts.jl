
using Pkg, DrWatson

@quickactivate "StatisticalRethinkingTuring"

indir = projectdir("notebooks")
outdir = projectdir("scripts")
!isdir(outdir) && mkdir(outdir)

function copy_file(
  fromfile::AbstractString, fromdir::AbstractString, todir::AbstractString
)

  nblines = readlines(joinpath(fromdir, fromfile))
  outfilename = joinpath(todir, fromfile)

  isfile(outfilename) && rm(outfilename)
  outfile = open(outfilename, "w")

  for i in nblines
    #println(i)
    if length(i) > 0 && i[1] == '#'
      continue
    else
      write(outfile, i)
      write(outfile, "\n")
    end
  end

  close(outfile)

end

cd(indir) do
  nbsubdirs = readdir()
  for nbsubdir in nbsubdirs
    if isdir(nbsubdir) && nbsubdir !== "scripts"
      if !isdir(joinpath(outdir, nbsubdir))
        mkdir(joinpath(outdir, nbsubdir))
      end
      # Find all notebooks in nbsubdir, skip subdir "intros"
      if !(nbsubdir == "intros")
        nbs = readdir(nbsubdir)
        println("$(nbsubdir): $(nbs)\n")
        # Copy the notebooks to the scripts dir
        for nb in nbs
          nb == ".DS_Store" && continue
          copy_file(nb, joinpath(indir, nbsubdir), joinpath(outdir, nbsubdir))
        end
      else # Handle subdirs in "intros"
        nb_intro_dirs = readdir(nbsubdir)
        println("$(nbsubdir): $(nb_intro_dirs)\n")
        # Copy the notebooks to the scripts dir
        for nb_intro_dir in nb_intro_dirs
          nb_intro_dir == ".DS_Store" && continue
          intro_nbs = readdir(joinpath(indir, nbsubdir, nb_intro_dir))
          println("$(nbsubdir)/$(nb_intro_dir): $(intro_nbs)\n")
          for nb in intro_nbs
            nb == ".DS_Store" && continue
            #println("Checking for $(joinpath(outdir, nbsubdir, nb_intro_dir))")
            !isdir(joinpath(outdir, nbsubdir, nb_intro_dir)) &&
              mkdir(joinpath(outdir, nbsubdir, nb_intro_dir))
            copy_file(nb, joinpath(indir, nbsubdir, nb_intro_dir),
              joinpath(outdir, nbsubdir, nb_intro_dir))
          end
        end
      end
    else
      println(("$(nbsubdir) ignored"))
    end
  end

end
