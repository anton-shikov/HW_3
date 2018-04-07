# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 10:24:30 2018

@author: anton
"""

from Bio import SeqIO
import collections
import numpy as np
import matplotlib.pyplot as plt

class Kmer_spectrum:
    k = 13
    sequence_name = ''
    kmer = ''
    
    def __init__(self, file_name, k, q):
        self.sequence_name = file_name
        self.kmer_quality = q
        self.kmer_size = k
        self.kmer_dict = collections.defaultdict(int)
    
    def increase(self, k_mer):
        self.kmer = k_mer
        self.kmer_dict[self.kmer] += 1
        
    def kmer_analyse(self):
           for record in SeqIO.parse(self.sequence_name, 'fastq'):
               seq = str(record.seq)
               seq_leng = len(record.seq)
               annotation_list = record.letter_annotations["phred_quality"]
               self.check_flag = 0
               for index in range (seq_leng-self.kmer_size+1):
                   for qual_check in range (index, index+self.kmer_size):
                       if annotation_list[qual_check] < self.kmer_quality:
                         self.check_flag = 1
                         break
                   if self.check_flag != 1:
                        current_kmer = seq[index:(index+self.kmer_size)]
                        self.increase(current_kmer)
                   self.check_flag = 0
           return(self.kmer_dict)
    
    def array_build (self, k_mer_dict):
        self.ready_dict = k_mer_dict
        self.unique_nums = list(set([self.ready_dict[key] for key in self.ready_dict.keys()]))
        self.counter_column = []
        for num in self.unique_nums:
            counter = 0
            for mer in self.ready_dict:
                if self.ready_dict[mer] == num:
                    counter += 1
            self.counter_column.append(counter)
        self.total_array = np.array(list(zip(self.unique_nums,self.counter_column)))
        return(self.total_array)
        
    def visualize_spec (self, k_mer_array, max_x, max_y):
        self.xmax, self.ymax = max_x, max_y
        self.array_to_viz = k_mer_array
        plt.bar(self.array_to_viz[:,0],  self.array_to_viz[:,1], align='center', color='#B00000', alpha=0.65)
        plt.xlabel('k-mer abundance')
        plt.ylabel('# of k-mers with that abundance')
        plt.title('k-mer spectrum: {}'.format(self.sequence_name))
        plt.axis([0, self.xmax, 0, self.ymax])
        plt.show()
        
    def genome_complex (self, noise_border, k_mer_array):
        self.border = noise_border
        self.raw_array = k_mer_array
        self.Array_to_analyse = np.delete(self.raw_array, [i for i in range (noise_border)], 0)
        summ = 0
        for k in range(len(self.Array_to_analyse)):
            summ += self.Array_to_analyse[k,0]*self.Array_to_analyse[k,1]
        return(self.summ/self.kmer_size)
        
Kmer = Kmer_spectrum ('C:\\Downloads\\test_kmer.fastq', 15, 20)
my_array = Kmer.array_build(Kmer.kmer_analyse())
Kmer.visualize_spec(my_array, 370, 40000)
Kmer.genome_complex(50, my_array)





