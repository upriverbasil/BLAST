import collections

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
	def __init__(self, start_pos, positions):
		self.start_pos = start_pos
		self.positions = positions
	def __eq__(self, other):
		if isinstance(other, Node):
			return self.start_pos == other.start_pos and self.positions == other.positions
		return False
	
	def __repr__(self):
		return f'Kmer start position:{self.start_pos} where it can be found:{self.positions}'

class BLAST:
	def __init__(self, database, query, word_length,HSSP = 1):
		self.database = database
		self.word_length = word_length
		self.kmers = {}
		self.hash_table = {'A':1,'G':2,'C':3,'T':4}
		self.kmers_hash_scores = {}
		self.binary_tree_array = []
		self.query = query
		self.HSSP = HSSP
		self.possible_strings = {}
		self.hit_kmers = []
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
			hash_score+=(self.hash_table[kmer[i]]*(4**i))
		return hash_score

	def kmers_positions(self):
		for i in range(0,len(self.database)-self.word_length+1):
			seq = self.database[i:i+self.word_length]
			if(seq in self.kmers):
				self.kmers[seq].append(i)
			else:
				self.kmers[seq] = [i]
		print(self.kmers)

	def kmers_hash_table(self):
		for i in self.kmers:
			self.kmers_hash_scores[i] = self.hash_code(i)
		print(self.kmers_hash_scores)

	def make_binary_tree(self):
		for i in self.kmers_hash_scores:
			self.binary_tree_array.append(Node(self.kmers_hash_scores[i], self.kmers[i]))
		self.binary_tree_array = sorted(self.binary_tree_array)
		print(self.binary_tree_array)

	def match_kmer_binary_tree(self):
		l = []
		for i in range(0,len(self.query)-self.word_length+1):
			seq = self.query[i:i+self.word_length]
			self.possible_strings = {seq}
			
			self.recur(seq,0,self.HSSP,0)
			print(self.possible_strings)
			for j in self.possible_strings:
				code = self.hash_code(j)
				idx = self.binary_search(code)
				if(idx!=None):
					l.append(K_Mer(i,idx.positions))
		self.hit_kmers = l
		print(self.hit_kmers)

if __name__ == '__main__':
	a = BLAST("ACTACT","ACG",3)
	a.kmers_positions()
	a.kmers_hash_table()
	a.make_binary_tree()
	print(a.binary_search(57))
	print(a.binary_search(78))
	print(a.binary_search(35))
	a.match_kmer_binary_tree()
