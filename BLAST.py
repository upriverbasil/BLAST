import collections
from os import stat
import numpy as np
from Bio import SeqIO
import math

class Node:
	def __init__(self, hash_score, positions):
		self.hash_score = hash_score
		self.positions = positions
	def __eq__(self, other):
		if isinstance(other, Node):
			return self.hash_score == other.hash_score
		return False
	def __cmp__(self,other):
		return self.hash_score-other.hash_score
	def __lt__(self, other):
		return self.hash_score<other.hash_score
	def __repr__(self):
		return f'Node hash_score:{self.hash_score} Position:{self.positions}'

class K_Mer:
	def __init__(self, start_pos, positions, sequence):
		self.start_pos = start_pos
		self.positions = positions
		self.sequence = sequence
	def __eq__(self, other):
		if isinstance(other, Node):
			return self.start_pos == other.start_pos and self.positions == other.positions and self.sequence == sequence
		return False

	def __repr__(self):
		return f'Kmer Sequence: {self.sequence} start position:{self.start_pos} where it can be found:{self.positions}'

class BLAST:
	def __init__(self, database, query, word_length,HSSP, insertion, deletion,mismatch,mscore):
		self.db_name = database
		self.insertion = insertion
		self.deletion = deletion
		self.mismatch = mismatch
		self.mscore = mscore
		self.database = ""
		self.word_length = word_length
		self.kmers = {}
		self.hash_table = {'A':0,'G':2,'C':1,'T':3}
		self.kmers_hash_scores = {}
		self.binary_tree_array = []
		self.query = query
		self.HSSP = HSSP
		self.possible_strings = {}
		self.hit_kmers = []

	def run(self):
		fasta_sequences = SeqIO.parse(open(self.db_name),'fasta')
		final_results_array = []
		for fasta in fasta_sequences:
			keys = (fasta.description).split(" ")
			key = keys[0]
			seq = str(fasta.seq)
			seq=seq.replace("N","")
			self.database = seq
			target = self.query
			self.kmers_positions()
			self.kmers_hash_table()
			self.make_binary_tree()
			self.match_kmer_binary_tree()
			print("q length -> ",len(target), " db_len -> ", len(self.database), " key-> ", key)
			key_res = self.smith_waterman(self.database,target,self.word_length,key)
			print(key_res, len(key_res))
			for i in key_res:
				final_results_array.append(key_res)
		final_results_array = np.array(final_results_array)
		print(final_results_array.shape, " <final shape")
		temp = final_results_array[final_results_array[:,0].argsort()]
		print(temp)
		final_results_array = temp

	def binary_search(self, val):
		low = 0
		high = len(self.binary_tree_array)-1
		while(low<=high):
			mid = int((low+high)/2)
			if(self.binary_tree_array[mid].hash_score==val):
				return self.binary_tree_array[mid]
			elif(self.binary_tree_array[mid].hash_score>val):
				high = mid-1
			else:
				low = mid+1
		return None

	def recur(self,s, mut, max, i):
		if(mut>=max or i>=len(s)):
			self.possible_strings.add(s)
		else:
			if(s[i]=="A"):
				self.recur(s,mut,max,i+1)
				self.recur(s[:i]+"G"+s[i+1:],mut+1,max,i+1)
				self.recur(s[:i]+"T"+s[i+1:],mut+1,max,i+1)
				self.recur(s[:i]+"C"+s[i+1:],mut+1,max,i+1)
			elif(s[i]=="G"):
				self.recur(s,mut,max,i+1)
				self.recur(s[:i]+"A"+s[i+1:],mut+1,max,i+1)
				self.recur(s[:i]+"T"+s[i+1:],mut+1,max,i+1)
				self.recur(s[:i]+"C"+s[i+1:],mut+1,max,i+1)
			elif(s[i]=="T"):
				self.recur(s,mut,max,i+1)
				self.recur(s[:i]+"A"+s[i+1:],mut+1,max,i+1)
				self.recur(s[:i]+"G"+s[i+1:],mut+1,max,i+1)
				self.recur(s[:i]+"C"+s[i+1:],mut+1,max,i+1)
			elif(s[i]=="C"):
				self.recur(s,mut,max,i+1)
				self.recur(s[:i]+"G"+s[i+1:],mut+1,max,i+1)
				self.recur(s[:i]+"T"+s[i+1:],mut+1,max,i+1)
				self.recur(s[:i]+"A"+s[i+1:],mut+1,max,i+1)

	def hash_code(self, kmer):
		hash_score = 0
		word_length = len(kmer)
		for i in range(0,word_length) :
			hash_score+=(self.hash_table[kmer[word_length-i-1]]*(4**i))

		return hash_score

	def kmers_positions(self):
		for i in range(0,len(self.database)-self.word_length+1):
			seq = self.database[i:i+self.word_length]
			if(seq in self.kmers):
				self.kmers[seq].append(i)
			else:
				self.kmers[seq] = [i]
		#print(self.kmers)

	def kmers_hash_table(self):
		for i in self.kmers:
			self.kmers_hash_scores[i] = self.hash_code(i)
		#print(self.kmers_hash_scores)

	def make_binary_tree(self):
		for i in self.kmers_hash_scores:
			self.binary_tree_array.append(Node(self.kmers_hash_scores[i], self.kmers[i]))
		self.binary_tree_array = sorted(self.binary_tree_array)
		#print(self.binary_tree_array)

	def match_kmer_binary_tree(self):
		l = []
		for i in range(0,len(self.query)-self.word_length+1):
			seq = self.query[i:i+self.word_length]
			self.possible_strings = {seq}

			self.recur(seq,0,self.HSSP,0)
			#print(self.possible_strings)
			for j in self.possible_strings:
				code = self.hash_code(j)
				idx = self.binary_search(code)
				#print(j,code)
				if(idx!=None):
					l.append(K_Mer(i,idx.positions,j))
		self.hit_kmers = l
		#print(self.hit_kmers)


	'''
	S - Alignment score
	K, lamba - Parameters
	m - length of database
	n - length of query sequence
	'''
	def statistics(self, n, S, m, lamba = 1.09861, K = 0.3333):
		bit_score = (lamba*S-math.log(K))/math.log(2)
		e_value = (m*n)/(2**bit_score)
		p_value = 2**(-1*bit_score)
		return [bit_score, e_value, p_value]
		# print("E-Value : "+str(e_value))
		# print("P-Value : "+str(p_value))
		# print("bit_score : "+str(bit_score))

	def smith_waterman(self,d,t,k, key):
		matrix= np.zeros((len(t)+1,len(d)+1))
		insertion= self.insertion
		deletion = self.deletion
		mismatch=self.deletion
		mscore=self.mscore
		#scoring the matrix based on the scores defined above
		for i in range(1,len(t)+1):
			for j in range(1,len(d)+1):
				matrix[i,j] = max(matrix[i][j-1] + insertion, matrix[i-1][j] + deletion, matrix[i-1][j-1] + (1 if t[i-1] == d[j-1] else +mismatch),0)

		pos=[]
		#extracting the previosuly generated alignment starting positions
		for i in self.hit_kmers:
			temp=[]
			temp.append(i.start_pos)
			for j in i.positions:
				temp.append(j)
			pos.append(temp)
		matches=[]
		start_positions = []
		quers=[]
		scores=[]
		#iterating over the matrix from each alignment starting position and tracing back the alignment
		print(len(pos))
		for i in pos:
			for j in range(0,len(i)-1):
				score=0
				match=""
				q=""
				x=i[0]+1
				y=i[j+1]+1
				iter=1
				cur=matrix[x][y]
				if t[x-1]==d[y-1]:
					match+=d[y-1]
					q+=t[x-1]
					score+=mscore
				else:
					match+=d[y-1]
					score+=insertion
					q+="-"
				#break condition 1 = we reach the end of the matrix
				while x< len(t)+1 and y<len(d):
					iter+=1
					if x!= len(t):
						right= matrix[x][y+1]
						bot= matrix[x+1][y]
						diag= matrix[x+1][y+1]
					else:
						right= matrix[x][y+1]
						bot= 0
						diag= 0
					cur=matrix[x][y]
					values=[right,bot,diag]
					direc = max(right,bot,diag)
					index_min = max(range(len(values)), key=values.__getitem__)
					#break condition 2 = length of alignment is divisible by k and further movement is not favorable i.e all indexes have lower score than present index in the matrix
					if (iter)%k==0:
						if cur<direc:
							continue
						else:
							break
					cur=direc
					#moving to the right
					if index_min==0:
						x=x
						y=y+1
						q+="-"
						match+=d[y-1]
						score+=insertion
					#moving to the bottom
					elif index_min==1:
						x=x+1
						y=y
						match+="-"
						q+=t[x-1]
						score+=deletion
					#moving diagnal
					else:
						x=x+1
						y=y+1
						match+=d[y-1]
						q+=t[x-1]
						score+=mscore
				min_thresh = max(1, 0.1*len(t))
				if(len(q)>min_thresh):
					matches.append(match)
					quers.append(q)
					start_positions.append(i[j+1])
					scores.append(score)
		#print("query",quers)
		#print(len(quers), len(scores), matrix.shape)
		#print(scores[450], quers[450], matches[450])
		#print("database",matches)
		#print(matrix)
		N = 10
		res = sorted(range(len(scores)), key = lambda sub: scores[sub])#[-N:] 
		res = res[::-1]
		enc_set = set()
		n_counter = 0
		final_ans = []
		for i in res:
			query_match = quers[i]
			db_match = matches[i]
			st_pos = start_positions[i]
			if st_pos not in enc_set:
				n_counter+=1
				if n_counter>N:
					break
				match_score = scores[i]
				stats = self.statistics(len(t), match_score, len(d))
				bit_score = stats[0]
				e_value = stats[1]
				p_value = stats[2]
				temp = [p_value, e_value, bit_score, st_pos, query_match, db_match, key]
				final_ans.append(temp)
				enc_set.add(st_pos)
		return final_ans
		#print(max(scores),max(range(len(scores)), key=scores.__getitem__))


if __name__ == '__main__':
	# fasta_sequences = SeqIO.parse(open("sequence.fasta"),'fasta')
	# seq=""
	# for fasta in fasta_sequences:
	# 	seq+=str(fasta.seq)
	# seq=seq.replace("N","")
	# #print(seq)
	# database= seq
	target = "GCCTATACAGTTGAACTCGGTACAGAAGTAAATGAGTTCGCCTGTGTTGTGGCAGATGCTGTCATAAAAACTTTGCAACCAGTATCTGAATTACTTACACCACTGGGCATTGATTTAGATGAGTGGAGTATGGCTACATACTACTTATTTGATGAGTCTGGTGAGTTTAAATTGGCTTC"
	# print(len(target))
	# k=3
	# a = BLAST(database,target,k)
	# a.kmers_positions()
	# a.kmers_hash_table()
	# a.make_binary_tree()
	# a.match_kmer_binary_tree()
	# a.smith_waterman(database,target,k)
	blast = BLAST("db_seq.fasta", target, 3, 1, -1,-1,-1,1)
	blast.run()
