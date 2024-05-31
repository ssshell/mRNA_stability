
import sys


class LoadAnnotation:
	"""load original CDS annotation file"""

	def __init__(self, filename_originalAnnotation):
		self.filename_originalAnnotation_ = filename_originalAnnotation
		self.gene_ = []
		self.originalStart_ = {}
		self.originalEnd_ = {}
		self.strand_ = {}

	def get_annotation(self):
		for line in self.filename_originalAnnotation_:
			line = line.rstrip().split("\t")
			gene = line[0]
			self.gene_.append(gene)
			self.originalStart_[gene] = int(line[1])
			self.originalEnd_[gene] = int(line[2])
			self.strand_[gene] = line[3]


class SplitAnnotation:
	"""get coordinates of partial CDS region.
	target region defined as 5'end, middle or 3'end, each is about 1/3 of entire CDS region."""

	def __init__(self, filename_geneLength, dict_original_start, dict_original_end, dict_strand, target_region):
		self.filename_geneLength_ = filename_geneLength
		self.gene_ = []
		self.gene_length_ = {}
		self.originalStart_ = dict_original_start
		self.originalEnd_ = dict_original_end
		self.strand_ = dict_strand
		self.partialStart_ = {}
		self.partialEnd_ = {}
		self.target_region_ = target_region

	def get_geneLength(self):
		next(self.filename_geneLength_)
		for line in self.filename_geneLength_:
			line = line.rstrip().split("\t")
			gene = line[0]
			self.gene_.append(gene)
			gene_length = int(line[1])
			self.gene_length_[gene] = gene_length

	def get_splitAnnotation(self):
		for gene in self.gene_:
			gene_length = self.gene_length_[gene]
			third_length = round(gene_length / 3)
			if self.strand_[gene] == "+":
				if self.target_region_ == "5p":
					partial_start = self.originalStart_[gene]
					partial_end = self.originalStart_[gene] + third_length - 1
					# print(gene, partial_start, third_length, partial_end)
				elif self.target_region_ == "mid":
					partial_start = self.originalStart_[gene] + third_length
					partial_end = self.originalStart_[gene] + third_length + third_length - 1
					# print(gene, partial_start, third_length, partial_end)
				else:
					partial_start = self.originalStart_[gene] + third_length + third_length
					partial_end = self.originalEnd_[gene]
					# print(gene, partial_start, third_length, partial_end)
			else:
				if self.target_region_ == "3p":
					partial_start = self.originalStart_[gene]
					partial_end = self.originalStart_[gene] + third_length - 1
					# print(gene, partial_start, third_length, partial_end)
				elif self.target_region_ == "mid":
					partial_start = self.originalStart_[gene] + third_length
					partial_end = self.originalStart_[gene] + third_length + third_length - 1
					# print(gene, partial_start, third_length, partial_end)
				else:
					partial_start = self.originalStart_[gene] + third_length + third_length
					partial_end = self.originalEnd_[gene]
					# print(gene, partial_start, third_length, partial_end)
			
			self.partialStart_[gene] = partial_start
			self.partialEnd_[gene] = partial_end
			# print(gene, partial_start, third_length, partial_end)


def main():
	f_originalAnnotation = open(sys.argv[1], 'r')
	original_annotation = LoadAnnotation(f_originalAnnotation)
	original_annotation.get_annotation()
	
	############## test session
	# print(original_annotation.gene_)
	# print(original_annotation.originalStart_)
	# print(original_annotation.originalEnd_)
	# print(original_annotation.strand_)
	############## end of test session

	####### get split region annotation (5'end, middle, 3'end)
	f_geneLength = open(sys.argv[2], 'r')
	split_target_region = sys.argv[3]
	split_annotation = SplitAnnotation(f_geneLength, original_annotation.originalStart_, original_annotation.originalEnd_, 
		original_annotation.strand_, split_target_region)
	split_annotation.get_geneLength()
	split_annotation.get_splitAnnotation()
	
	for gene in split_annotation.gene_:
		print(gene + "\t" + str(split_annotation.partialStart_[gene]) + "\t" + str(split_annotation.partialEnd_[gene]) + \
			"\t" + original_annotation.strand_[gene])

	############## test session
	# print(split_annotation.gene_length_)
	############## end of test session
	

if __name__ == "__main__":
	main()
