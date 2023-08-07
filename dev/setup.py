import yaml
import os


app_dir = os.getcwd()
if 'dev' not in app_dir:
    if 'test_data' not in app_dir:
        app_dir = os.path.join(app_dir, 'dev')
    else:
        app_dir = os.path.join(os.path.dirname(app_dir), 'dev')

config_file = open(os.path.join(app_dir, "config.yml"))
CONFIG = yaml.safe_load(config_file)
config_file.close()

vocab_file = open(os.path.join(app_dir, "vocabulary.yml"))
VOCABULARY = yaml.safe_load(vocab_file)
vocab_file.close()
