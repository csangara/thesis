import pandas as pd
import seaborn as sns
import numpy as np
import os
import matplotlib.pyplot as plt

##### DOWNSAMPLING SYNTHVISIUM #####
def plot_heatmap_synthvisium(model_text=False):
    methods = ['cell2location', 'music', 'RCTD', 'spotlight', 'stereoscope']
    path = r"D:\Work (Yr 2 Sem 1)\Thesis\downsampling\brain_cortex\rep1"
    
    fig, axs = plt.subplots(1, 5, figsize=(21,7))
    plt.subplots_adjust(wspace=0.25)
    
    for k in range(5):
        method = methods[k]
        
        # These methods sometimes fail, also has qacct file info
        if method in ['music', 'RCTD', 'spotlight']:
            qacct_filepath =  os.path.join(path, "logs")
            qacct_file = [file for file in os.listdir(qacct_filepath) \
                          if (file.endswith(".qacct") and file.startswith(method))]
                
            # Getting information from qacct file
            filename = os.path.join(qacct_filepath, qacct_file[0])
            info_dict = {i:{} for i in range(1,17)}
            foi = ["exit_status", "ru_wallclock", "maxvmem"] # fields of interest
            with open(filename, "r") as f:
                for line in f:
                    line = line.strip()
                    if line.startswith("="):
                        continue
                    title, content = line.split(maxsplit=1)
                    if title == "taskid":
                        taskid = int(content)
                    
                    elif title in foi:
                        info_dict[taskid][title] = content
            
            info_df = pd.DataFrame(info_dict).transpose()
            info_df['time'] = [int(x.rstrip("s"))/60 for x in info_df['ru_wallclock'].values] # entire jobtime
            mask = (info_df['exit_status'] == "1").values.reshape((4,4)) # tells whether task was successfule
        else:
            # for c2l and stereoscope, all tasks were successful
            mask=None
            
        # sns.heatmap(info_df['time'].values.reshape((4, 4)), mask=mask, center=20, annot=True)
        file_dir = os.path.join(path, "results", methods[k])
        
        files = [file for file in os.listdir(file_dir) if file.endswith(".txt") ]
        times = {i:0 for i in range(1,17)}
        model_time = 0
        for file in files:
            if method in ['music', 'RCTD', 'spotlight']:
                """
                with open(os.path.join(file_dir, file)) as f:
                    lines = f.readlines()
                    time = float(lines[-1].split()[1])
                    n = int(file.split("_")[2])
                    times[n] = time
                """
                times = info_df['time']
            else:
                with open(os.path.join(file_dir, file)) as f:
                    if 'model' not in file:
                        time = float(f.readline().strip())
                        n = int(file.split("_")[2])
                        times[n] = time/60
                    else:
                        model_time = float(f.readline().strip())/60
                        
        times_series = pd.Series(times).values.reshape((4,4))
        sns.heatmap(times_series, mask=mask, annot=True, square=True, linewidth=0.01,
                    xticklabels=[1000,5000,10000,15000], yticklabels=[100,1000,5000,10000],
                    fmt='.3g', ax=axs[k], cbar_kws={'label': 'Time (min)', 'orientation':'horizontal'},
                    cmap='hot_r')
        axs[k].set(xlabel="Genes", ylabel="Spots", title=methods[k])
        
        if model_text:
            model_time_text = "Model: {:.2f} min".format(model_time) if model_time else ""
            axs[k].text(2, 6.5, model_time_text, ha="center")
    
    extra_text = "_model" if model_text else ""
    plt.savefig(os.path.join(path, "plots\\all_times") + extra_text + ".png", bbox_inches='tight')

# plot_heatmap_synthvisium(True)

##### DOWNSAMPLING SCREF #####
def plot_heatmap_scref(model_text=False):
    methods = ['cell2location', 'music', 'RCTD', 'spotlight', 'stereoscope']
    path = r"D:\Work (Yr 2 Sem 1)\Thesis\downsampling\brain_cortex\rep1"
    
    fig, axs = plt.subplots(1, 5, figsize=(21,7))
    plt.subplots_adjust(wspace=0.25)
    
    for k in range(5):
        method = methods[k]
        
        # These methods sometimes fail, also has qacct file info
        if method in ['music', 'RCTD', 'spotlight']:
            qacct_filepath =  os.path.join(path, "logs_sc")
            qacct_file = [file for file in os.listdir(qacct_filepath) \
                          if (file.endswith(".qacct") and file.startswith(method))]
                
            # Getting information from qacct file
            filename = os.path.join(qacct_filepath, qacct_file[0])
            info_dict = {i:{} for i in range(1, 4)}
            foi = ["exit_status", "ru_wallclock", "maxvmem"] # fields of interest
            with open(filename, "r") as f:
                for line in f:
                    line = line.strip()
                    if line.startswith("="):
                        continue
                    title, content = line.split(maxsplit=1)
                    if title == "taskid":
                        taskid = int(content)
                    
                    elif title in foi:
                        info_dict[taskid][title] = content
            
            info_df = pd.DataFrame(info_dict).transpose()
            info_df['time'] = [int(x.rstrip("s"))/60 for x in info_df['ru_wallclock'].values] # entire jobtime
            mask = (info_df['exit_status'] == "1").values.reshape((1,3)) # tells whether task was successfule
            
            
        else:
            # for c2l and stereoscope, all tasks were successful
            mask=None
        
        # sns.heatmap(info_df['time'].values.reshape((4, 4)), mask=mask, center=20, annot=True)
        file_dir = os.path.join(path, "results_sc", methods[k])
        
        files = [file for file in os.listdir(file_dir) if file.endswith(".txt") ]
        times = {i:0 for i in range(1,4)}
        fit_times = np.array([])
        for file in files:
            if method in ['music', 'RCTD', 'spotlight']:
                # use job time, not just deconv
                times = info_df['time']
            else:
                with open(os.path.join(file_dir, file)) as f:
                    if 'model' not in file:
                        time = float(f.readline().strip())
                        n = int(file.split("_")[2])
                        times[n] += time/60
                        fit_times = np.append(fit_times, time)
                    else:
                        model_time = float(f.readline().strip())
                        n = int(file.split("_")[1])
                        times[n] += model_time/60
                        
        times_series = pd.Series(times).values.reshape((1,3))
        sns.heatmap(times_series, mask=mask, annot=True, square=True, linewidth=0.01,
                    xticklabels=[1000, 5000, 10000], yticklabels=['All'],
                    fmt='.3g', ax=axs[k], cbar_kws={'label': 'Time (min)', 'orientation':'horizontal'},
                    cmap='hot_r')
        axs[k].set(xlabel="Cells", ylabel="Genes", title=methods[k])
        
        if model_text:
            fit_time_text = "Fitting time: {:.1f} ± {:.1f} min".format(\
                            np.mean(fit_times)/60, np.std(fit_times)/60) if fit_times.size > 0 else ""
            axs[k].text(1.5, 2.8, fit_time_text, ha="center")
            
    extra_text = "_modelfitting" if model_text else ""
    plt.savefig(os.path.join(path, "plots\\all_times_sc") + extra_text + ".png", bbox_inches='tight')

plot_heatmap_scref(True)