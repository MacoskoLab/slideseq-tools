import subprocess
import re
import sys
from time import sleep, time
from math import ceil
import os
import json

# TODO if continue and ignore folder, will still register |ls|>1 so won't ever die.

# from new_submit_to_taskrunner import call_to_taskrunner
# for i in $(seq 50); do echo '["qsub", "-l", "h_vmem=1m", "-pe", "smp", "1", "-binding", "linear:1", "-l", "h_rt=00:30:00", "-N", "sleep", "/broad/macosko/jlanglie/05_sysadmin/just_sleep.sh"]$$' > $RANDOM; done
def subproc_out2arr(subproc_out):
    return list(filter(lambda x: x!="",
                            subproc_out.decode().split("\n")))

def subproc_res(cmd_array):
    process = subprocess.Popen(cmd_array, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # can't do wait, need communicate so no deadlock https://stackoverflow.com/questions/1445627/how-can-i-find-out-why-subprocess-popen-wait-waits-forever-if-stdout-pipe
    out, err = process.communicate()
    return subproc_out2arr(out), subproc_out2arr(err), process.returncode

def generate_rand(random_chars=12, alphabet="0123456789abcdefghijklmnopqrstuvwxyz"):
    r = random.SystemRandom()
    return ''.join([r.choice(alphabet) for i in range(random_chars)])

def main():
    MAX_NUM_JOBS = 950

    all_submitted_jobs = []
    already_submitted_jobs = set()

    num_times_nothing_to_do = 0

    if len(sys.argv) != 3:
        print("Please provide two arguments: script path and output_dir!")
        sys.exit()
        
    scriptpath = sys.argv[1]
    output_dir = sys.argv[2]
    taskrunner_dir = re.sub("/$", "", output_dir)+"/tmp_taskrunner/"

    while True:

        # put at top so if continue in the loop, still sleeps not a fast while loop
        # TODO make 30
        sleep(45)  # TODOO
        print("->")
        # oldest one on top
        out_ls, err_ls, status = subproc_res(["ls", "-rt1", taskrunner_dir])

        if err_ls != [] or status!=0:
            print("ERROR IN LS: ")
            print(err_ls)
            print(status)
            print(out_ls)
            continue

        out_ls = list(filter(lambda x: x!="done", out_ls))



        out_stat, err_stat, status = subproc_res(["qstat"])
        if err_stat != [] or status!=0:
            # just skip to get past, can do the set intersection below but seems complicated
            # ['critical error: getgrgid(1015) failed: Numerical result out of range']
            print("ERROR IN qstat: ... continuing")
            print(err_stat)
            print(status)
            print(out_stat)
            continue

        # # sometimes weird oneoff qstat errors, do twice and get unique
        # out_stat_2, err_stat, status = subproc_res(["qstat"])
        # if err_stat != [] or status!= 0 :
        #     import pdb; pdb.set_trace()
        #     print("ERROR IN qstat: ")
        #     print(err_ls)
        #     print(status)
        #     print(err_stat)
        #     continue

        # # only take lines in both, should be the same
        # out_stat = list(set(out_stat) & set(out_stat_2))

        # remove headers
        out_stat = list(filter(lambda l: len(l)>0 and not re.match("^(job-ID)|(---)", l), out_stat))

        # len(l)>0 above so [0] is safe
        running_jnumbers = [ln.strip().split()[0] for ln in out_stat]
        running_jnumbers = list(filter(lambda l: l!="", running_jnumbers))


        curr_num_jobs = len(out_stat)



        if len(out_ls) == 0:
            print("No files=no new jobs")

            # if there are no qstats running or all qstats aren't in jobs I submitted, wait 15 minutes then kill
            if len(set(running_jnumbers) & set(all_submitted_jobs)) == 0: # TODOOO PUT BACK
            #if len(set(running_jnumbers)) == 0: # TODOOO PUT BACK
                print("AND NO JOBS WE'RE WAITING ON *{}".format(num_times_nothing_to_do))
                num_times_nothing_to_do+=1
            if num_times_nothing_to_do > 25: #90: # 45 mintes/30 second intervals
                print("Welp it was a good run, but enough is enough (also there's nothin' else to do) :-)")
                sys.exit(0)
        else:

            num_times_nothing_to_do = 0
            print("HAS {} FILES/JOBS TO DO:".format(len(out_ls)))
            print(out_ls)


            if curr_num_jobs < MAX_NUM_JOBS:
                # worry is if multiple concurrent pipelines both try to add aggressibely
                # If nowhere near max_num_jobs, then do up to 1/5 of remaining space at one time
                # up to 75 at one time to be safe, maybe
                # if near max_num_jobs, then up to 10 at a time

                # designed this so worst case, if have 4-5 concurrent jobs,
                # should never go over. If have 6+ where all hit it in a
                # couple of cycles with 1k jobs, before qstat can catch up, hard to say
                # should catch on submission error !=0
                # https://www.wolframalpha.com/input/?i=plot+x%2B6*Piecewise%5B%7B%7B%28min%2875%2C+%281000-x%29%2F5%29%29%2C++1000-x+%3E+320%7D%2C+%7B%28min%2810%2C+%281000-x%29%2F5%29%29%2C+1000-x%3C%3D320%7D%7D%5D+for+0%3C%3Dx%3C%3D1000
                max_at_once = 70# 90
                if MAX_NUM_JOBS - curr_num_jobs < 200:
                    max_at_once = 20

                # TODO what if goes after 900 and then is negative
                # caughr above? probably
                add_up_to = max(0, min(max_at_once, ceil((MAX_NUM_JOBS-curr_num_jobs)/5)))

                num_adding = min(len(out_ls), add_up_to)

                print("Has space for {} new jobs ({} running), ... adding up to {} out of {}".format(MAX_NUM_JOBS-curr_num_jobs, curr_num_jobs, add_up_to, len(out_ls)))
                
                num_run = 0
                for fn in out_ls:
                    if num_run >= num_adding:
                        break
                    
                    print("Running file {}".format(fn))

                    if os.path.isdir(taskrunner_dir+fn):  
                        print("SKIPPING DIR {}".format(fn))
                        continue

                    if fn in already_submitted_jobs:
                        print("ALREADY SUBMITTED! SKIPPING (maybe mv needs to propagate)")
                        continue

                    with open(taskrunner_dir+fn) as f:
                        lines = f.readlines()
                        # short circuit logic means this will be ok on [0]
                        if len(lines) < 1  or re.search("\$\$$", lines[0]) is None:
                            print("FILE {} is empty/not finished with @@@@!. Waiting for later!".format(fn))
                            continue

                        if len(lines) > 1:
                            print("FILE {} is has 2 lines!!...writing to new temp files I guess".format(fn))
                            continue
                            # if try to recover could hit an infinite loop where makes more and more, shouldn't happen anyway
                            # try to recover if somehow have collision by writing new lines to random file
                            # for l in lines[1:]
                            #     with open(taskrunner_dir+"00_MORE_THAN_2__"+int(time())+generate_rand(), "w") as newf:
                            #         f.write(l)
                            #         f.close()

                        this_qsub = lines[0] 

                        # remove ending $$ as EOL
                        this_qsub = re.sub("\$\$$", "", this_qsub.strip())

                        """
                        # split by @@@
                        # this_qsub_arr = this_qsub.strip().split("@@@")
                        """
                        this_qsub_arr = json.loads(this_qsub.strip())

                        out_qsub, err_qsub, status = subproc_res(this_qsub_arr)

                        # will also catch if can't submit bc overfull so can retry but don't want to depend too much

                        if status!=0:
                            print("ERROR RUNNING JOB {}, retrying later".format(fn))
                            print(out_qsub)
                            print(status)
                            print(err_qsub)
                            continue


                        this_jid = ""
                        out_qsub = (out_qsub[0] if len(out_qsub) >=1 else "")
                        if "Your job" in out_qsub: 
                            out_qsub_arr = out_qsub.split()
                            # Your job 21710621 ("sleep") has been submitted
                            this_jid = out_qsub_arr[2] if len(this_qsub_arr) >=3 else "INVALID"
                            all_submitted_jobs.append(this_jid)
                        print("RUNNING JID {}".format(this_jid))
                        # import pdb; pdb.set_trace()

                        # move to done
                        if not os.path.exists(taskrunner_dir+"done/"):
                            os.makedirs(taskrunner_dir +"done/")
                        os.rename(taskrunner_dir+fn, taskrunner_dir+"done/"+fn+"__"+this_jid+"__"+str(int(time())))

                        already_submitted_jobs.add(fn)
                        num_run+=1
            else:
                print("FULL, let's wait")
            print("<-")
            # import pdb; pdb.set_trace()


# todo
# keep list of job ids and parse qstat. If all qstats are not in submission list ("" not in !) and has been 1/2 an hour, kill
# if files waiting, sleep 15-30. If no files, sleep 90 seconds
# add try catch so never stops on weird error



    


# add name for this exact one, if name doesn't exist

if __name__ == "__main__":
    main()

