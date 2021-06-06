import pandas as pd
from PIL import Image, ImageDraw

#### INITIALIZATION #####
path = r"D:\Work (Yr 2 Sem 1)\Thesis"

dataset_types = ["artificial_uniform_distinct", "artificial_diverse_distinct", "artificial_uniform_overlap", "artificial_diverse_overlap",
                  "artificial_dominant_celltype_diverse", "artificial_partially_dominant_celltype_diverse",
                  "artificial_dominant_rare_celltype_diverse", "artificial_regional_rare_celltype_diverse"]
datasets = ['brain_cortex_generation', 'cerebellum_cell_generation', 'cerebellum_nucleus_generation',
              'hippocampus_generation', 'kidney_generation', 'pbmc_generation', 'scc_p5_generation']
ref=True
ref_text = "_ref" if ref else ""

# Each plot has a different layout :(
dims_dict = {"accuracy": [111, 368, 87, 257, 15, 14],
             "sensitivity":[111, 368, 87, 257, 15, 15],
             "specificity":[111, 367, 87, 256.5, 16.5, 15],
             "precision":[98, 357, 87, 256, 14.5, 16],
             "F1": [97.5, 357, 88, 256, 14.5, 16.5],
             "prc": [111, 367, 88, 257, 16.5, 16]}
if ref:
    dims_dict["precision"] = [111, 367, 87, 256, 16.5, 15]
    
colors_dict = {"c2l":"#f8766d", "MuSiC":"#a3a500", "RCTD":"#00bf7d",\
               "SPOT":"#00b0f6", "stereo":"#e76bf3", "tie":"#a1a1a1"}

##### MAIN CODE #####
for metric in dims_dict:
    # Read output from R
    df = pd.read_csv(path+ r'\Misc\best_values\best_values_' + metric + '.tsv', sep="\t")
    
    # Create new section to count ties
    df['comb'] = df['dataset'] + "," + df['dataset_type']
    reps = df['comb'].value_counts()
    
    # Store best performing method (or tie) in matrix
    results = [[0]*8 for i in range(7)]
    for i, dataset in enumerate(datasets):
        for j, dt in enumerate(dataset_types):
            if reps[dataset+","+dt] > 1:
                results[i][j] = "tie"
            else:
                results[i][j] = df[(df['dataset']==dataset) & \
                                   (df['dataset_type']==dt)]["method"].values[0]
    
    # Calculate facet width and height
    d = dims_dict[metric]
    start_x, start_y = d[0], d[2]
    height, width = d[3]-d[2], d[1]-d[0]
    space_x, space_y = d[4], d[5]
    
    # Draw rectangles on image
    im = Image.open(path + r'\plots\metrics_facet\all_' + metric + '_facet' + ref_text + '.png')
    im = im.convert("RGB")
    img1 = ImageDraw.Draw(im)  
    for i in range(7):
        for j in range(8):
            # Some math
            x = j*(width+space_x)+start_x
            y = i*(height+space_y)+start_y
            shape = [(x, y), (x+width, y+height)]
            
            # Outline is the best performing-method
            img1.rectangle(shape, outline=colors_dict[results[i][j]], width=3)
    im.save(path + r'\plots\metrics_facet\all_' + metric + '_facet' + ref_text + '_boxed.png')
    
#### EXTRA: count best method ####
perf_list = []
for metric in dims_dict:
    print(metric)
    # Read output from R
    df = pd.read_csv(path+ r'\Misc\best_values\best_values_' + metric + '.tsv', sep="\t")
    
    # Create new section to count ties
    df['comb'] = df['dataset'] + "," + df['dataset_type']
    reps = df['comb'].value_counts()
    
    # Store best performing method (or tie) in matrix
    methods_dict = {method:0 for method in colors_dict.keys()}
    for i, dataset in enumerate(datasets):
        for j, dt in enumerate(dataset_types):
            if reps[dataset+","+dt] > 1:
                methods_dict["tie"] += 1
            else:
                best_method = df[(df['dataset']==dataset) & \
                                   (df['dataset_type']==dt)]["method"].values[0]
                methods_dict[best_method] += 1
    print(methods_dict)
    perf_list.append(pd.Series(methods_dict, name=metric))

pd.concat(perf_list, axis=1).to_csv(path + r'\Misc\best_values\best_values_count.csv')