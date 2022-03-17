
using Random

include("usefulFunctions.jl");

intervalExperiment = nothing;
# intervalExperiment = (1, 4000);

fileNames = datasetNames();
for fileName in fileNames
    for strategy in 1:4
        fid = open(string("script_", fileName,"_strategy_",strategy),"w");
        write(fid,"#!/bin/bash\n");
        write(fid,"#SBATCH --job-name=", fileName, "_strategy_", string(strategy), "\n");
        if intervalExperiment!=nothing
            write(fid,"#SBATCH --array=", string(intervalExperiment[1]), "-", string(intervalExperiment[2]), "\n");
        end;

        write(fid,"#SBATCH --time=4-04:00:00\n");
        write(fid,"#SBATCH -p shared\n");
        write(fid,"#SBATCH --qos shared\n");
        write(fid,"#SBATCH -C cl2s\n");

        write(fid,"\n");
        write(fid,"module load julia/1.4.2\n");
        write(fid,"julia experimentsDoME.jl ", fileName, " ", string(strategy), " \$SLURM_ARRAY_TASK_ID\n");
        close(fid);
    end;
end;

fid = open("executeAll","w");
for numFileName in 1:length(fileNames)
    fileName = fileNames[numFileName];
    for strategy in 1:4
        write(fid,"sbatch script_", fileName, "_strategy_", string(strategy), "\n");
    end;
end;
close(fid);

fid = open("cancelAll","w");
for fileName in fileNames
    for strategy in 1:4
        write(fid,"scancel -n ", fileName, "_strategy_", string(strategy), "\n");
    end;
end;
close(fid);
