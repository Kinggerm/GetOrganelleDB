#! /usr/bin/env python
import os
import sys
import subprocess
from multiprocessing import Pool
from copy import deepcopy
from optparse import OptionParser
from GetOrganelleLib.sam_parser import MapRecords, mapping_gap_info_from_coverage_dict
from GetOrganelleLib.seq_parser import SequenceList
import random
import time
# import intervals as itv
import shutil
import math


def get_options(print_title=""):
    parser = OptionParser(usage="generate_bowtie2_seed.py", description=print_title)
    parser.add_option("-o", dest="output",
                      help="(Required) output directory.")
    parser.add_option("-d", dest="input_dir",
                      help="Input directory that contains a list of fasta files (.fasta).")
    parser.add_option("-l", dest="largest_gap", type=int, default=2000,
                      help="Threshold of largest gap. Default: %default")
    parser.add_option("-k", dest="key_samples",
                      help="Use-defined key samples (split by commas or each sample per line in a file) to accept.")
    parser.add_option("-P", dest="prior_samples",
                      help="Use-defined prior samples (split by commas or each sample per line in a file) to check, "
                           "without guarantee of acceptance.")
    parser.add_option("-E", dest="exclude_samples",
                      help="Use-defined excluded samples (split by commas or each sample per line in a file), "
                           "might be included in the query but not database. ")
    parser.add_option("-N", dest="num_of_initial_searching", default=200, type=int,
                      help="Files to initial search (all-by-all bowtie2) for each rounds. Suggested < 500. "
                           "Default: %default")
    parser.add_option("-n", dest="num_of_expanding_searching", default=500, type=int,
                      help="Files to expanding search (all-by-all bowtie2) for each rounds. Default: %default")
    parser.add_option("-R", dest="random_start_rounds", default=10, type=int,
                      help="Rounds of searching for initial files, "
                           "after which files remained would be the final survivors. Default: %default")
    parser.add_option("-t", "--threads", dest="threads", default=2, type=int,
                      help="Threads for bowtie2. Default: %default")
    parser.add_option("--read-len", dest="simulating_read_len", default=100, type=int,
                      help="Read length for simulating reads. Default: %default")
    parser.add_option("--read-jump", dest="simulating_read_jump", default=2, type=int,
                      help="Jumping step for simulating reads. Default: %default")
    parser.add_option("--not-circular", dest="circular", default=True, action="store_false",
                      help="Sequences are not circular.")
    parser.add_option("--random-seed", dest="random_seed", default=12345, type=int,
                      help="Random seed.")
    parser.add_option("--verbose", dest="verbose_log", default=False, action="store_true",
                      help="Verbose logging.")
    parser.add_option("--restart", dest="restart", default=False, action="store_true",
                      help="Force restart if files existed. Default: %default")
    parser.add_option("--keep-temp", dest="keep_temp", default=False, action="store_true",
                      help="Keep temp files.")
    options, argv = parser.parse_args()
    if not ((len(argv) or options.input_dir) and options.output):
        sys.stdout.write("Insufficient arguments:\n")
        parser.print_help()
        sys.exit()
    return options, argv


def making_bowtie2_index(fasta_file, seed_index_base, random_seed=12345, log_handler=None, verbose_log=False):
    if os.path.exists(seed_index_base + ".rev.1.bt2l"):
        pass
    else:
        if log_handler:
            log_handler.info("Making seed - bowtie2 index ...")
        build_seed_index = subprocess.Popen("bowtie2-build --seed " + str(random_seed) + " --large-index " +
                                            fasta_file + " " + seed_index_base,
                                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        if log_handler and verbose_log:
            log_handler.info("bowtie2-build --seed " + str(random_seed) + " --large-index " +
                             fasta_file + " " + seed_index_base)
        output, err = build_seed_index.communicate()
        if "unrecognized option" in output.decode("utf8"):
            build_seed_index = subprocess.Popen("bowtie2-build --seed " + str(random_seed) + " " +
                                                fasta_file + " " + seed_index_base,
                                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            if log_handler and verbose_log:
                log_handler.info("bowtie2-build --seed " + str(random_seed) + " " + fasta_file + " " + seed_index_base)
            output, err = build_seed_index.communicate()
        if "(ERR)" in output.decode("utf8") or "Error:" in output.decode("utf8"):
            raise Exception(output.decode("utf8"))
        if log_handler and verbose_log:
            log_handler.info("\n" + output.decode("utf8").strip())
        if log_handler:
            log_handler.info("Making seed - bowtie2 index finished.")


def generate_fq_file(from_fasta_file, out_base, is_circular, sim_read_len, sim_read_jump_size):
    from_fq_file = os.path.join(out_base, os.path.basename(from_fasta_file[:-5] + "fq"))
    from_record = SequenceList(from_fasta_file).sequences[0]
    from_sequence = from_record.seq
    if set(from_sequence).issubset({"N", "-", "n"}):
        return 0
    from_seq_len = len(from_sequence)
    if not os.path.exists(from_fq_file):
        with open(from_fq_file + ".Temp", "w") as output_handler:
            if is_circular:
                from_sequence += from_sequence[:sim_read_len - 1]
            for go_base in range(0, len(from_sequence) - sim_read_len + 1, sim_read_jump_size):
                output_handler.write("".join(["@", str(go_base), "\n",
                                              from_sequence[go_base: go_base + sim_read_len], "\n",
                                              "+", str(go_base), "\n",
                                              "G" * sim_read_len, "\n"]))
        os.rename(from_fq_file + ".Temp", from_fq_file)
    return from_seq_len


def get_largest_gap(from_fasta_file, to_fasta_file, from_seq_len, out_base, cal_database_gap_from_map_info=False,
                    sim_read_len=100, sim_read_jump_size=2, is_circular=True,
                    keep_temp=False, bowtie2_other_options="", verbose_log=False,
                    log_handler=None, threads=1, random_seed=12345):
    from_fq_file = os.path.join(out_base, os.path.basename(from_fasta_file[:-5] + "fq"))
    seed_index_base = os.path.join(out_base, os.path.basename(to_fasta_file))
    # making_bowtie2_index(to_fasta_file, seed_index_base)
    # mapping
    if log_handler:
        log_handler.info("Mapping reads to seed - bowtie2 index ...")
    mapped_fq_f = os.path.join(out_base, os.path.basename(from_fq_file).replace(".fq", "") + "--to--" +
                               os.path.basename(to_fasta_file).replace(".fasta", "") + ".fq")
    mapped_sam_f = os.path.join(out_base, os.path.basename(from_fq_file).replace(".fq", "") + "--to--" +
                                os.path.basename(to_fasta_file).replace(".fasta", "") + ".sam")

    def do_bowtie2():
        this_command = "bowtie2 --seed " + str(random_seed) + " --mm -p " + str(threads) + " " \
                       + bowtie2_other_options + \
                       " -x " + seed_index_base + " -U " + from_fq_file + " -S " + mapped_sam_f + ".Temp" + \
                       " --no-unal --omit-sec-seq" + " --al " + mapped_fq_f + ".Temp"
        make_seed_bowtie2 = subprocess.Popen(this_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        if verbose_log and log_handler:
            log_handler.info(this_command)
        output, err = make_seed_bowtie2.communicate()
        if "(ERR)" in output.decode("utf8") or "Error:" in output.decode("utf8"):
            raise Exception(output.decode("utf8"))
        if verbose_log and log_handler:
            log_handler.info("\n" + output.decode("utf8").strip())
        if os.path.exists(mapped_fq_f + ".Temp"):
            os.rename(mapped_sam_f + ".Temp", mapped_sam_f)
            os.rename(mapped_fq_f + ".Temp", mapped_fq_f)
            if log_handler:
                log_handler.info("Mapping finished.")
        else:
            raise Exception("Cannot find bowtie2 result!")
    if not (os.path.exists(mapped_sam_f) and os.path.exists(mapped_fq_f)):
        do_bowtie2()
        # this_command = "bowtie2 --seed " + str(random_seed) + " --mm -p " + str(threads) + " " \
        #                + bowtie2_other_options + \
        #                " -x " + seed_index_base + " -U " + from_fq_file + " -S " + mapped_sam_f + \
        #                " --no-unal --omit-sec-seq" + " --al " + mapped_fq_f
        # make_seed_bowtie2 = subprocess.Popen(this_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        # if verbose_log and log_handler:
        #     log_handler.info(this_command)
        # output, err = make_seed_bowtie2.communicate()
        # if "(ERR)" in output.decode("utf8") or "Error:" in output.decode("utf8"):
        #     raise Exception(output.decode("utf8"))
        # if verbose_log and log_handler:
        #     log_handler.info("\n" + output.decode("utf8").strip())
        # if os.path.exists(mapped_sam_f):
        #     if log_handler:
        #         log_handler.info("Mapping finished.")
        # else:
        #     raise Exception("Cannot find bowtie2 result!")
    else:
        if log_handler:
            log_handler.info("Mapping reads to seed - bowtie2 index ... skipped")
    if cal_database_gap_from_map_info:
        # could somehow save time
        try:
            mapping_records = MapRecords(mapped_sam_f, log_handler=log_handler)
            mapping_records.update_coverages()
            if mapping_records.coverages:
                database_gaps = mapping_gap_info_from_coverage_dict(mapping_records.coverages)[2]
                largest_db_gap = database_gaps[0] if database_gaps else 0
            else:
                largest_db_gap = list(mapping_records.references.values())[0]["real_len"]
        except (ValueError, IndexError) as e:
            print("Regenerate sam file: ", mapped_sam_f)
            do_bowtie2()
            mapping_records = MapRecords(mapped_sam_f, log_handler=log_handler)
            mapping_records.update_coverages()
            if mapping_records.coverages:
                database_gaps = mapping_gap_info_from_coverage_dict(mapping_records.coverages)[2]
                largest_db_gap = database_gaps[0] if database_gaps else 0
            else:
                largest_db_gap = list(mapping_records.references.values())[0]["real_len"]
            # raise e
    else:
        largest_db_gap = None
    seq_ids = []
    with open(mapped_fq_f) as open_fq_handler:
        for line in open_fq_handler:
            seq_ids.append(int(line[1:].strip()))
            for foo in range(3):
                open_fq_handler.readline()
    seq_ids.sort()
    gap_start = False
    gap_end = False
    gaps = []
    if not seq_ids:
        gaps.append(from_seq_len)
    else:
        go_to_id = 0
        if seq_ids[go_to_id]:
            gap_start = True
            gaps.append(seq_ids[go_to_id])
        while go_to_id < len(seq_ids):
            go_to_id += 1
            while go_to_id < len(seq_ids) and seq_ids[go_to_id] - seq_ids[go_to_id - 1] - sim_read_len<= 0:
                go_to_id += 1
            if go_to_id < len(seq_ids):
                gaps.append(seq_ids[go_to_id] - seq_ids[go_to_id - 1] - sim_read_len)
            else:
                if from_seq_len - seq_ids[go_to_id - 1] - sim_read_len > 0:
                    gap_end = True
                    gaps.append(from_seq_len - seq_ids[go_to_id - 1] - sim_read_len)
                else:
                    if gap_start:
                        gaps[0] = max(0, gaps[0] + from_seq_len - seq_ids[go_to_id - 1] - sim_read_len)
                        if not gaps[0]:
                            del gaps[0]
                break
        if gap_start and gap_end and is_circular:
            gaps[0] += gaps.pop(-1)
    if not keep_temp:
        os.remove(mapped_sam_f)
        os.remove(mapped_fq_f)
    ref_gap = sorted(gaps, reverse=True)[0] if gaps else 0
    return ref_gap, largest_db_gap


def pool_multiprocessing(target, iter_args, constant_args, num_process):
    try:
        # parse args
        if iter_args:
            if type(iter_args) in {list, tuple}:
                if type(iter_args[0]) not in {list, tuple}:
                    iter_args = [[each_arg] for each_arg in iter_args]
            else:
                sys.stderr.write("iter_args must be list/tuple!\n")
        else:
            return
        if constant_args:
            if type(constant_args) not in {list, tuple}:
                constant_args = [constant_args]
        else:
            constant_args = []
        results = []
        pool = Pool(processes=num_process)
        for this_arg in iter_args:
            results.append(pool.apply_async(target, tuple(list(this_arg) + list(constant_args))))
        pool.close()
        try:
            pool.join()
        except KeyboardInterrupt:
            pool.terminate()
            raise KeyboardInterrupt
        else:
            return [func_res.get() for func_res in results]
    except KeyboardInterrupt:
        raise KeyboardInterrupt


def select_backbone_fasta_files(output_dir, candidates, calling_network, threshold,
                                cal_database_gap_from_map_info=False,
                                user_designed_set=None,
                                sim_read_len=100, sim_read_jump_size=2, is_circular=True,
                                keep_temp=False, verbose_log=False,
                                log_handler=None, threads=1, random_seed=12345):
    """
    :param output_dir:
    :param candidates:
    :param calling_network:
    :param threshold:
    :param cal_database_gap_from_map_info: True only when all database files contain only one sequence
    :param user_designed_set:
    :param sim_read_len:
    :param sim_read_jump_size:
    :param is_circular:
    :param keep_temp:
    :param verbose_log:
    :param log_handler:
    :param threads:
    :param random_seed:
    :return:
    """
    time0 = time.time()
    discarded_samples = set()
    if not user_designed_set:
        user_designed_set = set()

    # generating fq and index
    query_seq_lengths = {}
    if threads == 1:
        for query_seq in candidates:
            if query_seq not in user_designed_set:
                query_seq_lengths[query_seq] = generate_fq_file(query_seq, output_dir, is_circular,
                                                                sim_read_len, sim_read_jump_size)
        for database_seq in user_designed_set:
            seed_index_base = os.path.join(output_dir, os.path.basename(database_seq))
            making_bowtie2_index(database_seq, seed_index_base, random_seed)
    else:
        iter_args = [query_seq for query_seq in candidates if query_seq not in user_designed_set]
        constant_args = [output_dir, is_circular, sim_read_len, sim_read_jump_size]
        seq_lens = pool_multiprocessing(generate_fq_file, iter_args, constant_args, threads)
        for go_q, query_seq in enumerate(iter_args):
            query_seq_lengths[query_seq] = seq_lens[go_q]
        iter_args = []
        for database_seq in user_designed_set:
            seed_index_base = os.path.join(output_dir, os.path.basename(database_seq))
            iter_args.append([database_seq, seed_index_base])
        pool_multiprocessing(making_bowtie2_index, iter_args, random_seed, threads)
    with open(os.path.join(output_dir, "invalid_seqs.txt"), "a") as invalid_seq_h:
        for query_seq in query_seq_lengths:
            if not query_seq_lengths[query_seq]:
                invalid_seq_h.write(query_seq + "\n")
                candidates.remove(query_seq)
    # map to user_designed set
    this_calling_network = {}

    def distribute_gap_to_network(here_iter_args, here_gaps, here_this_calling_network, here_calling_network):
        for go_res, (q_seq, db_seq, q_len) in enumerate(here_iter_args):
            if db_seq not in here_calling_network:
                here_calling_network[db_seq] = {}
                here_this_calling_network[db_seq] = {}
            elif db_seq not in here_this_calling_network:
                here_this_calling_network[db_seq] = {}
            if q_seq not in here_calling_network[db_seq]:
                here_ref_gap, here_db_gap = here_gaps[go_res]
                here_this_calling_network[db_seq][q_seq] = \
                    here_calling_network[db_seq][q_seq] = here_ref_gap
                if cal_database_gap_from_map_info:
                    if q_seq not in here_calling_network:
                        here_calling_network[q_seq] = {}
                    here_calling_network[q_seq][db_seq] = here_db_gap
            elif q_seq not in here_this_calling_network[db_seq]:
                here_this_calling_network[db_seq][q_seq] = here_calling_network[db_seq][q_seq]

    if user_designed_set:
        # all-to-designed
        if threads == 1:
            for query_seq in candidates:
                if query_seq not in user_designed_set:
                    for database_seq in user_designed_set:
                        if database_seq not in calling_network:
                            calling_network[database_seq] = {}
                            this_calling_network[database_seq] = {}
                        elif database_seq not in this_calling_network:
                            this_calling_network[database_seq] = {}
                        if query_seq not in calling_network[database_seq]:
                            ref_gap, db_gap = \
                                get_largest_gap(query_seq, database_seq, query_seq_lengths[query_seq], output_dir,
                                                cal_database_gap_from_map_info=cal_database_gap_from_map_info,
                                                sim_read_len=sim_read_len, sim_read_jump_size=sim_read_jump_size,
                                                is_circular=is_circular,
                                                keep_temp=keep_temp, verbose_log=verbose_log,
                                                log_handler=log_handler, threads=1, random_seed=random_seed)
                            this_calling_network[database_seq][query_seq] = \
                                calling_network[database_seq][query_seq] = ref_gap
                            if cal_database_gap_from_map_info:
                                if query_seq not in calling_network:
                                    calling_network[query_seq] = {}
                                calling_network[query_seq][database_seq] = db_gap
                        elif query_seq not in this_calling_network:
                            this_calling_network[database_seq][query_seq] = calling_network[database_seq][query_seq]
        else:
            iter_args = []
            for query_seq in candidates:
                if query_seq not in user_designed_set:
                    for database_seq in user_designed_set:
                        iter_args.append((query_seq, database_seq, query_seq_lengths[query_seq]))
            constant_args = [output_dir, cal_database_gap_from_map_info, sim_read_len, sim_read_jump_size,
                             is_circular, keep_temp, "", verbose_log, log_handler, 1, random_seed]
            gaps = pool_multiprocessing(get_largest_gap, iter_args, constant_args, threads)
            distribute_gap_to_network(iter_args, gaps, this_calling_network, calling_network)
        time1 = time.time()
        print("    candidates-to-designed cost: " + "%.4f" % (time1 - time0))
        # clear up
        for database_seq in this_calling_network:
            for query_seq in list(this_calling_network[database_seq]):
                if query_seq in this_calling_network[database_seq]:
                    if this_calling_network[database_seq][query_seq] > threshold:
                        del this_calling_network[database_seq][query_seq]
        # print(this_calling_network)
        # print(user_designed_set)
        for user_f in user_designed_set:
            for query_seq_reached in this_calling_network[user_f]:
                if query_seq_reached in this_calling_network:
                    print("query deleted: " + query_seq_reached)
                    del this_calling_network[query_seq_reached]
                if query_seq_reached in candidates:
                    candidates.remove(query_seq_reached)
                    discarded_samples.add(query_seq_reached)
        # print(candidates)
        time2 = time.time()
        print("    clearing up cost: " + "%.4f" % (time2 - time1))
    time2 = time.time()
    # survivors-to-survivors
    # making bowtie2 database
    s_t_s_cal_database_gap = True
    for query_seq in user_designed_set:
        candidates.remove(query_seq)
    if len(candidates) > 1:
        sorted_candidates = sorted(candidates)
        if threads == 1:
            for database_seq in sorted_candidates:
                seed_index_base = os.path.join(output_dir, os.path.basename(database_seq))
                making_bowtie2_index(database_seq, seed_index_base, random_seed)
            for query_seq in sorted_candidates:
                for database_seq in sorted_candidates:
                    if query_seq != database_seq:
                        if database_seq not in calling_network:
                            calling_network[database_seq] = {}
                            this_calling_network[database_seq] = {}
                        elif database_seq not in this_calling_network:
                            this_calling_network[database_seq] = {}
                        if query_seq not in calling_network[database_seq]:
                            ref_gap, db_gap = \
                                get_largest_gap(query_seq, database_seq, query_seq_lengths[query_seq], output_dir,
                                                cal_database_gap_from_map_info=s_t_s_cal_database_gap,
                                                sim_read_len=sim_read_len, sim_read_jump_size=sim_read_jump_size,
                                                is_circular=is_circular,
                                                keep_temp=keep_temp, verbose_log=verbose_log,
                                                log_handler=log_handler, threads=threads, random_seed=random_seed)
                            this_calling_network[database_seq][query_seq] = \
                                calling_network[database_seq][query_seq] = ref_gap
                            if s_t_s_cal_database_gap:
                                if query_seq not in calling_network:
                                    calling_network[query_seq] = {}
                                calling_network[query_seq][database_seq] = db_gap
                        elif query_seq not in this_calling_network:
                            this_calling_network[database_seq][query_seq] = calling_network[database_seq][query_seq]
        else:
            iter_args = []
            for database_seq in sorted_candidates:
                seed_index_base = os.path.join(output_dir, os.path.basename(database_seq))
                iter_args.append([database_seq, seed_index_base])
            pool_multiprocessing(making_bowtie2_index, iter_args, random_seed, threads)
            iter_args = []
            existed_pairs = []
            for query_seq in sorted_candidates:
                for database_seq in sorted_candidates:
                    if query_seq != database_seq:
                        if database_seq not in calling_network or query_seq not in calling_network[database_seq]:
                            iter_args.append((query_seq, database_seq, query_seq_lengths[query_seq]))
                        else:
                            existed_pairs.append((query_seq, database_seq, query_seq_lengths[query_seq]))
            constant_args = [output_dir, s_t_s_cal_database_gap, sim_read_len, sim_read_jump_size,
                             is_circular, keep_temp, "", verbose_log, log_handler, 1, random_seed]
            gaps = pool_multiprocessing(get_largest_gap, iter_args, constant_args, threads)
            distribute_gap_to_network(iter_args, gaps, this_calling_network, calling_network)
            distribute_gap_to_network(existed_pairs, [], this_calling_network, calling_network)
        time3 = time.time()
        print("    survivors-to-survivors(" + str(len(candidates)) + ") cost: " + "%.4f" % (time3 - time2))

        # clear up
        for database_seq in this_calling_network:
            # if database_seq == "fungus_mitogenome-20190507/Fungi__380714936__JQ346808.fasta":
            #     print("checking fungus_mitogenome-20190507/Fungi__380714936__JQ346808.fasta")
            #     print(this_calling_network[database_seq])
            # elif database_seq == "fungus_mitogenome-20190507/Fungi__380714957__JQ346809.fasta":
            #     print("checking fungus_mitogenome-20190507/Fungi__380714957__JQ346809.fasta")
            #     print(this_calling_network[database_seq])
            # elif database_seq == "fungus_mitogenome-20190507/Fungi__381211955__NC_017016.fasta":
            #     print("checking fungus_mitogenome-20190507/Fungi__381211955__NC_017016.fasta")
            #     print(this_calling_network[database_seq])
            for query_seq in list(this_calling_network[database_seq]):
                if query_seq in this_calling_network[database_seq]:
                    if this_calling_network[database_seq][query_seq] > threshold:
                        del this_calling_network[database_seq][query_seq]
        for user_f in user_designed_set:
            for query_seq_reached in this_calling_network[user_f]:
                if query_seq_reached in this_calling_network:
                    del this_calling_network[query_seq_reached]
                    discarded_samples.add(query_seq_reached)
        survivors = sorted(this_calling_network,
                           key=lambda x:
                           (-len(this_calling_network[x]),
                            *sorted([this_calling_network[x][q] for q in this_calling_network[x]], reverse=True),
                            x))
        # greedy search
        for check_sv in survivors:
            # if check_sv == "fungus_mitogenome-20190507/Fungi__380714936__JQ346808.fasta":
            #     print("checking fungus_mitogenome-20190507/Fungi__380714936__JQ346808.fasta")
            # elif check_sv == "fungus_mitogenome-20190507/Fungi__380714957__JQ346809.fasta":
            #     print("checking fungus_mitogenome-20190507/Fungi__380714957__JQ346809.fasta")
            # elif check_sv == "fungus_mitogenome-20190507/Fungi__381211955__NC_017016.fasta":
            #     print("checking fungus_mitogenome-20190507/Fungi__381211955__NC_017016.fasta")
            if check_sv not in user_designed_set and check_sv in this_calling_network:
                # if check_sv == "fungus_mitogenome-20190507/Fungi__380714936__JQ346808.fasta":
                #     print("checking fungus_mitogenome-20190507/Fungi__380714936__JQ346808.fasta")
                #     print(this_calling_network[check_sv])
                # elif check_sv == "fungus_mitogenome-20190507/Fungi__380714957__JQ346809.fasta":
                #     print("checking fungus_mitogenome-20190507/Fungi__380714957__JQ346809.fasta")
                #     print(this_calling_network[check_sv])
                # elif check_sv == "fungus_mitogenome-20190507/Fungi__381211955__NC_017016.fasta":
                #     print("checking fungus_mitogenome-20190507/Fungi__381211955__NC_017016.fasta")
                #     print(this_calling_network[check_sv])
                for query_seq_reached in this_calling_network[check_sv]:
                    # if check_sv == "fungus_mitogenome-20190507/Fungi__380714936__JQ346808.fasta":
                    #     if query_seq_reached == "fungus_mitogenome-20190507/Fungi__380714957__JQ346809.fasta":
                    #         print("809 in 808, and 809 in network ", query_seq_reached in this_calling_network)
                    # if check_sv == "fungus_mitogenome-20190507/Fungi__380714936__JQ346809.fasta":
                    #     if query_seq_reached == "fungus_mitogenome-20190507/Fungi__380714957__JQ346809.fasta":
                    #         print("808 in 809, and 808 in network ", query_seq_reached in this_calling_network)
                    # if check_sv == "checking fungus_mitogenome-20190507/Fungi__381211955__NC_017016.fasta":
                    #     if query_seq_reached == "fungus_mitogenome-20190507/Fungi__380714936__JQ346808.fasta":
                    #         print("808 in 016, and 808 in network ", query_seq_reached in this_calling_network)
                    if query_seq_reached in this_calling_network:
                        del this_calling_network[query_seq_reached]
                        discarded_samples.add(query_seq_reached)
        print("    discarded samples: " + str(len(discarded_samples)) + "; "
              "survivors: " + str(len(this_calling_network) - len(user_designed_set)))
    elif len(candidates) == 1:
        # if only one candidates left
        this_calling_network[candidates.pop()] = {}
        time3 = time.time()
        print("    survivors-to-survivors(" + str(len(candidates)) + ") cost: " + "%.4f" % (time3 - time2))
    else:
        # if no candidates
        time3 = time.time()
        print("    survivors-to-survivors(" + str(len(candidates)) + ") cost: " + "%.4f" % (time3 - time2))
    # exhausted search
    # connected parts
    # backup calling_network
    with open(os.path.join(output_dir, "calling_network.Temp"), "w") as writing_net:
        writing_net.write(str(calling_network))
    os.rename(os.path.join(output_dir, "calling_network.Temp"), os.path.join(output_dir, "calling_network.txt"))
    time4 = time.time()
    print("    clearing up cost: " + "%.4f" % (time4 - time3))
    return set(this_calling_network), discarded_samples


def generate_samples(p_samples, untested_samples, num_samples):
    if p_samples:
        samples_to_return = random.sample(sorted(p_samples), min(num_samples, len(p_samples)))
        for inside_check_s in samples_to_return:
            untested_samples.remove(inside_check_s)
            p_samples.remove(inside_check_s)
        if not p_samples and untested_samples and len(samples_to_return) < num_samples:
            extra_s = random.sample(sorted(untested_samples),
                                    min(num_samples - len(samples_to_return), len(untested_samples)))
            for inside_check_s in extra_s:
                untested_samples.remove(inside_check_s)
            samples_to_return.extend(extra_s)
    else:
        samples_to_return = random.sample(sorted(untested_samples), min(num_samples, len(untested_samples)))
        for inside_check_s in samples_to_return:
            untested_samples.remove(inside_check_s)
    return samples_to_return


def main():
    time_start = time.time()
    options, argv = get_options()
    try:
        random.seed(options.random_seed)
        if options.input_dir:
            all_samples = argv + [os.path.join(options.input_dir, fs)
                                  for fs in os.listdir(options.input_dir) if fs.endswith(".fasta")]
        else:
            all_samples = argv
        remaining_untested = set(all_samples)
        if options.key_samples:
            if os.path.isfile(options.key_samples):
                if options.key_samples.endswith(".fasta"):
                    user_designed_samples = {str(options.key_samples)}
                else:
                    user_designed_samples = set([user_s.strip()
                                                 for user_s in open(options.key_samples) if user_s.strip()])
            else:
                user_designed_samples = set(options.key_samples.split(","))
        else:
            user_designed_samples = set()
        if options.prior_samples:
            if os.path.isfile(options.prior_samples):
                prior_samples = set([user_s.strip() for user_s in open(options.prior_samples) if user_s.strip()])
            else:
                prior_samples = set(options.prior_samples.split(","))
        else:
            prior_samples = set()
        for check_p_s in prior_samples:
            assert check_p_s in remaining_untested, "Prior sample " + check_p_s + "  not found!"
        if options.restart and os.path.exists(options.output):
            shutil.rmtree(options.output)
        if not os.path.exists(options.output):
            os.mkdir(options.output)
        # load previous network if existed
        if os.path.exists(os.path.join(options.output, "calling_network.txt")):
            calling_network = eval(open(os.path.join(options.output, "calling_network.txt")).read())
        else:
            calling_network = {}
        # remove user-defined from untested

        for user_s in user_designed_samples:
            remaining_untested.remove(user_s)
            prior_samples.discard(user_s)
        time_preparing = time.time()
        print("Preparing cost: " + "%.4f" % (time_preparing - time_start))

        # initial searching
        initial_prior_samples = deepcopy(prior_samples)
        initial_remaining_untested = deepcopy(remaining_untested)
        initial_existed = [int(str(i_file.split("R")[1]).split(".")[0]) for i_file in os.listdir(options.output)
                           if i_file.startswith("initial_in_samples_R") and i_file.endswith(".txt")]
        if initial_existed:
            largest_initial = max(initial_existed)
            print(str(largest_initial) + " initial rounds existed!")
            initial_samples = set()
            with open(os.path.join(options.output, "initial_in_samples_R" + str(largest_initial) + ".txt")) as ini_h:
                for initial_s in ini_h:
                    if initial_s.strip():
                        initial_samples.add(initial_s.strip())
                        if initial_s.strip() not in user_designed_samples:
                            initial_remaining_untested.remove(initial_s.strip())
                        initial_prior_samples.discard(initial_s.strip())
            for e_f_n in initial_existed:
                with open(os.path.join(options.output, "initial_ex_samples_R" + str(e_f_n) + ".txt")) as dump_handler:
                    for discard_s in dump_handler:
                        if discard_s.strip():
                            initial_remaining_untested.remove(discard_s.strip())
                            initial_prior_samples.discard(discard_s.strip())
            initial_pseudo_prior_samples = deepcopy(prior_samples)
            initial_pseudo_remaining_untested = deepcopy(remaining_untested)
        else:
            largest_initial = 0
            initial_samples = set()
            initial_pseudo_prior_samples = set()
            initial_pseudo_remaining_untested = set()

        time_reading_initial = time.time()
        print("Reading initial existed cost: " + "%.4f" % (time_reading_initial - time_preparing))
        print("Initial samples: " + str(len(initial_samples)) + "; Untested: " + str(len(initial_remaining_untested)))
        for initial_count in range(options.random_start_rounds):
            round_start_time = time.time()
            if not initial_remaining_untested:
                break
            # to reproduce, make random.sample used.
            if initial_count < largest_initial:
                generate_samples(initial_pseudo_prior_samples, initial_pseudo_remaining_untested,
                                 options.num_of_initial_searching)
                continue
            samples_this_round = generate_samples(initial_prior_samples, initial_remaining_untested,
                                                  options.num_of_initial_searching)
            for in_check in samples_this_round:
                initial_samples.add(in_check)
            for in_check in user_designed_samples:
                initial_samples.add(in_check)
            initial_samples, discarded_this_round = select_backbone_fasta_files(
                options.output, deepcopy(initial_samples), calling_network,
                threshold=options.largest_gap,
                cal_database_gap_from_map_info=True,
                user_designed_set=user_designed_samples,
                sim_read_len=options.simulating_read_len, sim_read_jump_size=options.simulating_read_jump,
                is_circular=options.circular, keep_temp=options.keep_temp, verbose_log=options.verbose_log,
                log_handler=None, threads=options.threads, random_seed=options.random_seed)
            with open(os.path.join(options.output, "initial_in_samples_R" + str(initial_count + 1) + ".temp"), "w") as a_h:
                a_h.writelines([a_sam + "\n" for a_sam in sorted(initial_samples)])
            with open(os.path.join(options.output, "initial_ex_samples_R" + str(initial_count + 1) + ".temp"), "w") as d_h:
                d_h.writelines([d_sam + "\n" for d_sam in sorted(discarded_this_round)])
            os.rename(os.path.join(options.output, "initial_ex_samples_R" + str(initial_count + 1) + ".temp"),
                      os.path.join(options.output, "initial_ex_samples_R" + str(initial_count + 1) + ".txt"))
            os.rename(os.path.join(options.output, "initial_in_samples_R" + str(initial_count + 1) + ".temp"),
                      os.path.join(options.output, "initial_in_samples_R" + str(initial_count + 1) + ".txt"))
            round_end_time = time.time()
            print("Initial round: " + str(initial_count + 1) + "/" + str(options.random_start_rounds) +
                  "  Accumulated samples: " + str(len(initial_samples)) +
                  "  Cost: " + "%.4f" % (round_end_time - round_start_time))

        # extending searching
        num_accumulated_sample = len(initial_samples)
        for initial_s in initial_samples:
            if initial_s not in user_designed_samples:
                remaining_untested.remove(initial_s)
            prior_samples.discard(initial_s)
        extending_existed = [int(str(i_file.split("R")[1]).split(".")[0]) for i_file in os.listdir(options.output)
                             if i_file.startswith("extending_in_samples_R") and i_file.endswith(".txt")]
        # concatenate fasta
        # if not os.path.exists(os.path.join(options.output, "extending_in_samples_R0.fasta")):
        with open(os.path.join(options.output, "extending_in_samples_R0.fasta.Temp"), "w") as extending_fasta:
            for initial_s in sorted(initial_samples):
                with open(initial_s) as ini_fs:
                    extending_fasta.write(ini_fs.read().strip() + "\n")
        os.rename(os.path.join(options.output, "extending_in_samples_R0.fasta.Temp"),
                  os.path.join(options.output, "extending_in_samples_R0.fasta"))
        if extending_existed:
            largest_extending = max(extending_existed)
            print(str(largest_extending) + " extending rounds existed!")
            for i_f_n in extending_existed:
                with open(os.path.join(options.output, "extending_in_samples_R" + str(i_f_n) + ".txt")) as ini_h:
                    for extend_s in ini_h:
                        if extend_s.strip():
                            remaining_untested.remove(extend_s.strip())
                            prior_samples.discard(extend_s.strip())
                            num_accumulated_sample += 1
            for e_f_n in extending_existed:
                with open(os.path.join(options.output, "extending_ex_samples_R" + str(e_f_n) + ".txt")) as dump_handler:
                    for discard_s in dump_handler:
                        if discard_s.strip():
                            remaining_untested.remove(discard_s.strip())
                            prior_samples.discard(discard_s.strip())
            extending_pseudo_prior_samples = deepcopy(prior_samples)
            extending_pseudo_remaining_untested = deepcopy(remaining_untested)
        else:
            largest_extending = 0
            extending_pseudo_prior_samples = set()
            extending_pseudo_remaining_untested = set()
        print("Reading extending existed cost: " + "%.4f" % (time_reading_initial - time_preparing))
        num_remaining_rounds = int(math.ceil(len(remaining_untested) / options.num_of_expanding_searching))
        for extending_count in range(largest_extending + num_remaining_rounds):
            round_start_time = time.time()
            # to reproduce, make random.sample used.
            if extending_count < largest_extending:
                generate_samples(extending_pseudo_prior_samples, extending_pseudo_remaining_untested,
                                 options.num_of_expanding_searching)
                continue
            samples_this_round = generate_samples(prior_samples, remaining_untested, options.num_of_expanding_searching)
            previous_accepted_sample = \
                os.path.join(options.output, "extending_in_samples_R" + str(extending_count) + ".fasta")
            # helps under rerun mode; no effect on new runs;
            if previous_accepted_sample in calling_network:
                del calling_network[previous_accepted_sample]
            check_samples_this_round = {previous_accepted_sample}
            for in_check in samples_this_round:
                check_samples_this_round.add(in_check)
            extending_samples_this_round, discarded_this_round = select_backbone_fasta_files(
                options.output, check_samples_this_round, calling_network,
                threshold=options.largest_gap,
                cal_database_gap_from_map_info=False,
                user_designed_set={previous_accepted_sample},
                sim_read_len=options.simulating_read_len, sim_read_jump_size=options.simulating_read_jump,
                is_circular=options.circular, keep_temp=options.keep_temp, verbose_log=options.verbose_log,
                log_handler=None, threads=options.threads, random_seed=options.random_seed)
            # concatenate fasta
            if len(discarded_this_round) == len(samples_this_round):
                # previous_accepted_sample not changed
                base_name = os.path.join(options.output, "extending_in_samples_R")
                for postfix in ("", ".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l"):
                    shutil.copy(base_name + str(extending_count) + ".fasta" + postfix,
                                base_name + str(extending_count + 1) + ".fasta" + postfix)
            else:
                with open(os.path.join(options.output, "extending_in_samples_R" + str(extending_count + 1) + ".fasta.Temp"), "w") as extending_fasta:
                    for extend_s in sorted(extending_samples_this_round):
                        with open(extend_s) as ini_fs:
                            extending_fasta.write(ini_fs.read().strip() + "\n")
                os.rename(os.path.join(options.output, "extending_in_samples_R" + str(extending_count + 1) + ".fasta.Temp"),
                          os.path.join(options.output, "extending_in_samples_R" + str(extending_count + 1) + ".fasta"))
            # write out sample names
            extending_samples_this_round.remove(previous_accepted_sample)
            num_accumulated_sample += len(extending_samples_this_round)
            with open(os.path.join(options.output, "extending_in_samples_R" + str(extending_count + 1) + ".temp"), "w") as a_h:
                a_h.writelines([a_sam + "\n" for a_sam in sorted(extending_samples_this_round)])
            with open(os.path.join(options.output, "extending_ex_samples_R" + str(extending_count + 1) + ".temp"), "w") as d_h:
                d_h.writelines([d_sam + "\n" for d_sam in sorted(discarded_this_round)])
            os.rename(os.path.join(options.output, "extending_ex_samples_R" + str(extending_count + 1) + ".temp"),
                      os.path.join(options.output, "extending_ex_samples_R" + str(extending_count + 1) + ".txt"))
            os.rename(os.path.join(options.output, "extending_in_samples_R" + str(extending_count + 1) + ".temp"),
                      os.path.join(options.output, "extending_in_samples_R" + str(extending_count + 1) + ".txt"))
            round_end_time = time.time()
            print("Extending round: " + str(extending_count + 1) + "/" + str(largest_extending + num_remaining_rounds) +
                  "  Accumulated samples: " + str(num_accumulated_sample) +
                  "  Cost: " + "%.4f" % (round_end_time - round_start_time))
    except KeyboardInterrupt:
        pass
    print("Total cost: " + "%.4f" % (time.time() - time_start))


if __name__ == '__main__':
    main()

