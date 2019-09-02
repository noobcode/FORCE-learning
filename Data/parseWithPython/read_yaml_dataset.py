import yaml
import yamlordereddictloader

from yaml import load, dump, safe_load
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


#with open("../traj1.yaml") as f:
#	print("2")
#	yaml_data = yaml.load(f, Loader=yamlordereddictloader.Loader)
#	print("3")

#print(yaml.dump(yaml_data))

with open("../traj1.yaml") as f:
	yaml_data = load(f, Loader=Loader)

print(yaml.dump(yaml_data))
print(yaml_data)
