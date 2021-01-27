import json
import re
import random
import datetime

def generate_rand(random_chars=10, alphabet="0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"):
    r = random.SystemRandom()
    return ''.join([r.choice(alphabet) for i in range(random_chars)])


def call_to_taskrunner(output_dir, call_args):
	# really really don't want a collision so put very accurate date + random number
	fn = datetime.datetime.now().strftime("%d_%m_%Y__%H_%M_%S_%f")+"__"+generate_rand()
	#remove endng / if exists on output to sanitize
	with open("{}/tmp_taskrunner/{}".format(re.sub("/$", "", output_dir),
						fn), "w+") as f:
		# add $$ at end so the taskrunner knows is totally done
		f.write(json.dumps(call_args)+"$$")
			
