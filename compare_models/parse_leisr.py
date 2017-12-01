from phyphy import *
import os


json_directory = "inference/"
jsons = [x for x in os.listdir(json_directory) if x.endswith("json")]

fit_output = "model_fits.csv"
fitstring = "dataset,model,logl,AICc\n"

for json in jsons:
    print json
    p = HyPhyParser(json_directory + json)
    
    output = json.replace(".json", ".csv")
    p.extract_csv(json_directory + output)
    
    model = p.extract_model_names()[0]
    
    dataset = json.split(".")[0]
    outmodelname = json.split(".")[2]

    thisfit =  json.split(".")[0] + "," + outmodelname + "," + str(p.extract_model_logl(model)) + "," + str(p.extract_model_aicc(model)) + "\n"
    fitstring += thisfit

fitstring.strip()
with open(fit_output, "w") as f:
    f.write(fitstring)
    