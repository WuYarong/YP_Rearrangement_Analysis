#!/usr/bin/python
# -*- coding: utf-8 -*-


from __future__ import print_function
import sys
import collections


def read_IS_blast(file, len_min=2000, prefix='.IS.m6.out'):
	stat_IS = collections.defaultdict(lambda: collections.defaultdict(lambda: 0))
	partial = collections.defaultdict(lambda: [])
	partial_merge = collections.defaultdict(lambda: [])
	partial_left = collections.defaultdict(lambda: [])
	partial_left_minlen = collections.defaultdict(lambda: [])
	partial_left_maxlen = collections.defaultdict(lambda: [])
	IS_len = dict()

	with open(file) as f:
		for line in f.readlines():
			line = line.replace('\n', '')
			tmp = line.split('\t')
			tag = (tmp[0].split(':'))[0].replace(prefix, '')
			IS_len[tag] = int(tmp[-1])

			if int(tmp[3]) >= int(tmp[-1]) * 0.9:
				stat_IS[tag]['intact'] += 1
			else:
				stat_IS[tag]['partial'] += 1
				partial[tag].append([int(tmp[6]), int(tmp[7])])


	for key in partial.keys():
		for i in range(0, len(partial[key])-1):
			for j in range(i, len(partial[key])):
				if max(partial[key][i]) > max(partial[key][j]):
					if min(partial[key][i]) - max(partial[key][j]) <= len_min and (partial[key][i][1] - partial[key][i][0] + partial[key][j][1] - partial[key][j][0] + 2 ) >= IS_len[key] * 0.9:
						stat_IS[key]['partial'] -= 1
						# print key, '\t', partial[key][i], '\t', partial[key][j]
						partial_merge[key].append(partial[key][i])
						partial_merge[key].append(partial[key][j])

				elif max(partial[key][i]) < max(partial[key][j]):
					if min(partial[key][j]) - max(partial[key][i]) <= len_min and (partial[key][i][1] - partial[key][i][0] + partial[key][j][1] - partial[key][j][0]) >= IS_len[key] * 0.9:
						stat_IS[key]['partial'] -= 1
						partial_merge[key].append(partial[key][i])
						partial_merge[key].append(partial[key][j])

		for item in partial[key]:
			flag = 0
			for item2 in partial_merge[key]:
				if item[0] == item2[0]:
					flag = 1
			if flag == 0:
				partial_left[key].append(item)

		for item in partial_left[key]:
			if item[1] - item[0] + 1 < IS_len[key] * 0.5:
				stat_IS[key]['partial'] -= 1
				partial_left_minlen[key].append(item)
			else:
				partial_left_maxlen[key].append(item)

	print('strain\tIS_total_num\tIS_intact_num\tIS_partial_total_num\tincomplete_num\tmerge_num\tIS_partial_loc\tpartial_incomplete_more-half-len\tpartial_merge\tpartial_incomplete_less-half-len')
	for key in stat_IS.keys():
		print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(key, stat_IS[key]['intact'] + stat_IS[key]['partial'], stat_IS[key]['intact'], stat_IS[key]['partial'], len(partial_left_maxlen[key]), len(partial_merge[key])/2, partial[key], partial_left_maxlen[key], partial_merge[key], partial_left_minlen[key]))


if __name__ == '__main__':
	if len(sys.argv) < 2:
		print('Usage: python {} <grep.merge.blast.m6.out>'.format(sys.argv[0]))
		sys.exit(-1)

	read_IS_blast(sys.argv[1])
