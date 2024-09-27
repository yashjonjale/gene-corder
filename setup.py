import json
import os
import sys



config = {}

org_list = ["kdjsj","skdm"]

config['supported_organisms'] = {}
for org in org_list:
    config['supported_organisms'][org] = {}

config['objects'] = {}

#writing to json file
with open('config.json', 'w') as f:
    json.dump(config, f, indent=4)
