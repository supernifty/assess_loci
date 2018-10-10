#!/usr/bin/env python
'''
  given two groups of samples, assess their ability to be correctly classified using MS loci
'''

import argparse
import logging
import sys

import cyvcf2
import intervaltree

REMOVE_CHR=True

def populate_intervals(panels, vdx, vcf, name, group, pass_only, loci):
  logging.info('processing %s...', vcf)
  considered = filtered = 0
  found = [0] * len(panels)
  vcf_in = cyvcf2.VCF(vcf)
  for variant in vcf_in:
    # is pass
    if pass_only and variant.FILTER is not None:
      filtered += 1
      continue
    # is indel?
    if len(variant.REF) != len(variant.ALT[0]):
      considered += 1
      # check all the panels
      for pdx, panel in enumerate(panels):
        if REMOVE_CHR and variant.CHROM.startswith('chr'):
          chr = variant.CHROM[3:]
        else:
          chr = variant.CHROM
        if chr in panel and len(panel[chr][variant.POS]) > 0:
          found[pdx] += 1
          for interval in panel[chr][variant.POS]:
            logging.debug(interval)
            interval[2][0].append(vdx)
            loci[pdx].add((chr, interval))

  logging.info('processing %s: considered %i filtered %i found %s', vcf, considered, filtered, ' '.join([str(x) for x in found]))

def main(vcfs, names, groups, panels, pass_only):
  logging.info('starting with %i vcfs and %i groups: %i in group 0 and %i group 1...', len(vcfs), len(groups), len([x for x in groups if x == 0]), len([x for x in groups if x == 1]))

  # build interval trees
  panel_intervals = []
  for bed in panels:
    logging.info('processing %s...', bed)
    intervals = {}
    total = 0
    for idx, line in enumerate(open(bed, 'r')):
      fields = line.strip('\n').split('\t')
      if len(fields) < 3:
        logging.warn('skipped line %i in %s', idx + 1, bed)
        continue
      chr, start, finish = fields[:3]
      if len(fields) > 3:
        annot = fields[3]
      else:
        annot = ''
      if REMOVE_CHR and chr.startswith('chr'):
        chr = chr[3:]
      if chr not in intervals:
        intervals[chr] = intervaltree.IntervalTree()
      if len(intervals[chr][int(start):int(finish)]) > 0:
        pass #logging.debug('overlap at %s %s %s', chr, start, finish)
      intervals[chr][int(start):int(finish)] = ([], annot) # list of matching samples
      total += int(finish) - int(start)
      if (idx + 1) % 100000 == 0:
        logging.debug('%i lines...', idx + 1)
    panel_intervals.append(intervals)
    logging.info('processing %s: %i bases', bed, total)
    
  # assign affected variants to panels
  loci = []
  for _ in panels:
    loci.append(set())
  if len(names) == 1:
    names = [names[0]] * len(vcfs)
  for vdx, (vcf, name, group) in enumerate(zip(vcfs, names, groups)):
    populate_intervals(panel_intervals, vdx, vcf, name, group, pass_only, loci)

  # accuracy of each locus
  sys.stdout.write('Panel\tChr\tStart\tEnd\tAnnot\tTP\tTN\tFP\tFN\tSpecificity\tSensitivity\tAccuracy\n')
  for pdx, panel in enumerate(panels): # each panel
    for locus in loci[pdx]:
      # how informative is this locus?
      tp = tn = fp = fn = 0
      for idx, group in enumerate(groups):
        if idx in locus[1][2][0]: # found in interval
          if group == 0: # no good
            fp += 1
          else: # good
            tp += 1
        else: # not found in interval
          if group == 0: # good
            tn += 1
          else: # no good
            fn += 1
      sensitivity = tp / (tp + fn)
      specificity = tn / (tn + fp)
      accuracy = (tp + tn) / (tp + tn + fp + fn)
      sys.stdout.write('{panel}\t{chr}\t{start}\t{end}\t{annot}\t{tp}\t{tn}\t{fp}\t{fn}\t{specificity:.2f}\t{sensitivity:.2f}\t{accuracy:.2f}\n'.format(
        panel=panel,
        chr=locus[0],
        start=locus[1][0],
        end=locus[1][1],
        annot=locus[1][2][1],
        tp=tp,
        tn=tn,
        fp=fp,
        fn=fn,
        sensitivity=sensitivity,
        specificity=specificity,
        accuracy=accuracy
      ))

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess classifiability')
  parser.add_argument('--vcfs', required=True, nargs='+', help='vcfs to analyse')
  parser.add_argument('--names', required=True, nargs='+', help='sample names')
  parser.add_argument('--groups', required=True, nargs='+', type=int, help='which group (0 or 1)')
  parser.add_argument('--panels', required=True, nargs='+', help='panels (bed file format) to assess')
  parser.add_argument('--filter_pass', action='store_true', help='only pass calls')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.vcfs, args.names, args.groups, args.panels, args.filter_pass)
