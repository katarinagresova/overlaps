#!/usr/bin/env python3
from operator import attrgetter
from collections import defaultdict
from bisect import bisect_left
from datetime import date
import argparse

# debug prints only sums of weights
debug = False


class Sequence:
    """ Sequence holds Regions and its total score and finish point.
    """
    def __init__(self, weight, finish, reg):
        """ Creates Sequence
        :param weight: total weight of Regions
        :param finish: finish point of Regions
        :param reg:    Regions
        """
        self.weight = weight
        self.finish = finish
        self.regs = reg

    def __add__(self, region):
        """ Adds Region to current Sequence.
        :param region: Region to add
        :return: Sequence with added Region
        """
        return Sequence(
                    self.weight + region.weight,
                    region.finish,
                    (*self.regs, region)
                )

    def __gt__(self, other):
        """ Compare Sequences by weight.
        :param other: Sequence to compare with
        :return: true, if current Sequence has bigger score, false otherwise
        """
        return self.weight > other.weight


class Region:
    """ Region is area in nucleic acid sequence represented by one line in gff file.
    """
    def __init__(self, line):
        """ Creates Region.
        :param line: line from gff file
        """
        self.line = line
        parts = line.split()
        self.rid = parts[0] + parts[2] + parts[6]
        self.start = int(parts[3])
        self.finish = int(parts[4])
        self.weight = int(parts[5])


class Resolver:
    """ Resolver encapsulate functionality for resolving overlaps.
    """

    def __init__(self, input_file, is_sorted):
        """ Creates Resolver.
        :param input_file: path to input file in gff format
        :param is_sorted:  flat, that input file is sorted by start
        """
        self.input_file = input_file
        self.is_sorted = is_sorted

    @staticmethod
    def get_partial(routes):
        """ Finds the best non overlapping regions.
        :param routes: regions to analyse
        :return: score of the best Sequence
        """
        sequences = defaultdict(Sequence)
        routes.sort(key=attrgetter('finish'))
        count = 0

        # Extract start and finish points
        start = [f.start for f in routes]
        finish = [f.finish for f in routes]

        # Initial conditions
        sequences[0] = Sequence(0, 0, set())
        sequences[1] = Sequence(routes[0].weight, routes[0].finish, (routes[0],))

        # Find the best sequences
        for i in range(2, len(routes) + 1):
            new_seq = sequences[bisect_left(finish, start[i - 1])] + routes[i - 1]
            sequences[i] = max(sequences[i - 1], new_seq)

        # Print the best sequences
        if not debug:
            for reg in sequences[len(routes)].regs:
                print(reg.line, end='')
        if debug:
            for reg in sequences[len(routes)].regs:
                count += 1


        return sequences[len(routes)].weight, count

    def read_input(self):
        """ Reads input file, parses it and creates Region object for each line.
        :return: list of Regions
        """

        with open(self.input_file, 'r') as file:

            data = defaultdict(list)
            for line in file:
                # Skip empty lines, meta-data and comments
                if len(line) < 2 or line[0] == '#':
                    continue

                # Make Region object of current line
                f = Region(line)
                data[f.rid].append(f)

        return data

    def resolve(self):
        """ Finds the best non overlapping sequences in input file.
        """

        if not debug:
            print('##gff-version 3')
            print('##date %s' % date.today())

        data = self.read_input()
        # Find best overlap for each distinct chromosome, strand and type
        for rid in data.keys():
            # Holder for sequences
            sequences = list()
            # Top score so far
            top_weight = 0
            # Most far end so far
            most_far = 0
            count = 0

            # If input file was not sorted, we will sort it.
            # This will allow as to divide file into blocks.
            if not self.is_sorted:
                data[rid].sort(key=attrgetter('start'))

            for region in data[rid]:
                if len(sequences) > 0 and region.start > most_far:
                    top_weight_tmp, tmp_count = self.get_partial(sequences)
                    top_weight += top_weight_tmp
                    count += tmp_count
                    sequences.clear()

                most_far = max(region.finish, most_far)
                sequences.append(region)

            top_weight_tmp, tmp_count = self.get_partial(sequences)
            top_weight += top_weight_tmp
            count += tmp_count

            # Print score statistics
            if debug:
                print("# Type '{0}' scored: {1:10}".format(rid, top_weight))
                print(count)


if __name__ == '__main__':
    # Parse command line arguments
    ap = argparse.ArgumentParser()
    ap.add_argument("input", nargs=1, help="Path to input file in gff format.")
    ap.add_argument("-s", required=False, action='store_true', help="Flag, that input file is sorted by start.")
    args = ap.parse_args()
    # Resolve overlaps
    Resolver(args.input[0], args.s).resolve()
