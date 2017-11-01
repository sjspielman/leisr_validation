from phyphy import *
import os


json_directory = "inference/"
jsons = [x for x in os.listdir(json_directory) if x.endswith("json")]


for json in jsons:
    
    p = HyPhyParser(json_directory + json)
    
    output = json.replace(".json", ".csv")
    p.extract_csv(json_directory + output)