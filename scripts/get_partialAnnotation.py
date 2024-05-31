
import sys


class LoadAnnotation:
	"""load CDS annotation"""

	def __init__(self, filename_CDSannotation):
		self.filename_CDSannotation_ = filename_CDSannotation
		self.CDS_ = []
		self.CDS_start_ = {}
		self.CDS_end_ = {}
		self.CDS_strand_ = {}

	def get_CDSannotation(self):
		for line in self.filename_CDSannotation_:
			line = line.rstrip().split("\t")
			CDS = line[0]
			self.CDS_.append(CDS)
			self.CDS_start_[CDS] = int(line[1])
			self.CDS_end_[CDS] = int(line[2])
			self.CDS_strand_[CDS] = line[3]


class GetPartialAnnotation:
	"""get annotation of smaller portions of complete annotation
	either by percentage of length or number of nucleotides"""

	def __init__(self, list_CDS, dict_CDS_start, dict_CDS_end, dict_CDS_strand):
		self.CDS_ = list_CDS
		self.CDS_start_ = dict_CDS_start
		self.CDS_end_ = dict_CDS_end
		self.CDS_strand_ = dict_CDS_strand
		self.CDS_partial_coordinates_byPortion_ = {}
		self.CDS_partial_coordinates_byNucleotide_from5p_ = {}
		self.CDS_partial_coordinates_byNucleotide_from3p_ = {}

	def get_partial_byPortion(self, n_partial_portion_):
		"""slice by number of parts at the end"""

		for CDS_i in self.CDS_:
			self.CDS_partial_coordinates_byPortion_[CDS_i] = []
			CDS_length_i = self.CDS_end_[CDS_i] - self.CDS_start_[CDS_i] + 1
			portion_length = round(CDS_length_i / n_partial_portion_)
			self.CDS_partial_coordinates_byPortion_[CDS_i].append(self.CDS_start_[CDS_i])
			for j in range(n_partial_portion_ + 1):
				if len(self.CDS_partial_coordinates_byPortion_[CDS_i]) < n_partial_portion_:
					coordinate_j = self.CDS_start_[CDS_i] + (j + 1) * portion_length - 1
					# print(coordinate_j)
					self.CDS_partial_coordinates_byPortion_[CDS_i].append(coordinate_j)
			if self.CDS_end_[CDS_i] not in self.CDS_partial_coordinates_byPortion_[CDS_i]:
				self.CDS_partial_coordinates_byPortion_[CDS_i].append(self.CDS_end_[CDS_i])
			# print(CDS_i, CDS_length_i, portion_length, n_parts, self.CDS_start_[CDS_i], self.CDS_end_[CDS_i])
			# print(self.CDS_partial_coordinates_byPortion_[CDS_i])
			
	def get_partial_byNucleotide_from5p(self, n_partial_nucleotide_):
		"""slice by number of nucleotides are from 5'end starting position"""

		for CDS_i in self.CDS_:
			self.CDS_partial_coordinates_byNucleotide_from5p_[CDS_i] = []
			if self.CDS_strand_[CDS_i] == "+":
				self.CDS_partial_coordinates_byNucleotide_from5p_[CDS_i].append(self.CDS_start_[CDS_i])
				partial_end_i = self.CDS_start_[CDS_i] + n_partial_nucleotide_ - 1
				if partial_end_i >= self.CDS_end_[CDS_i]:
					self.CDS_partial_coordinates_byNucleotide_from5p_[CDS_i].append(self.CDS_end_[CDS_i])
				else:
					self.CDS_partial_coordinates_byNucleotide_from5p_[CDS_i].append(partial_end_i)
			else:
				partial_end_i = self.CDS_end_[CDS_i] - n_partial_nucleotide_ + 1
				if partial_end_i <= self.CDS_start_[CDS_i]:
					self.CDS_partial_coordinates_byNucleotide_from5p_[CDS_i].append(self.CDS_start_[CDS_i])
				else:
					self.CDS_partial_coordinates_byNucleotide_from5p_[CDS_i].append(partial_end_i)
				self.CDS_partial_coordinates_byNucleotide_from5p_[CDS_i].append(self.CDS_end_[CDS_i])
			# print(CDS_i, n_partial_nucleotide_, self.CDS_start_[CDS_i], self.CDS_end_[CDS_i])
			# print(self.CDS_partial_coordinates_byNucleotide_from5p_[CDS_i])

	def get_partial_byNucleotide_from3p(self, n_partial_nucleotide_):
		"""slice by number of nucleotides are from 3'end ending position"""

		for CDS_i in self.CDS_:
			self.CDS_partial_coordinates_byNucleotide_from3p_[CDS_i] = []
			if self.CDS_strand_[CDS_i] == "-":
				self.CDS_partial_coordinates_byNucleotide_from3p_[CDS_i].append(self.CDS_start_[CDS_i])
				partial_end_i = self.CDS_start_[CDS_i] + n_partial_nucleotide_ - 1
				if partial_end_i >= self.CDS_end_[CDS_i]:
					self.CDS_partial_coordinates_byNucleotide_from3p_[CDS_i].append(self.CDS_end_[CDS_i])
				else:
					self.CDS_partial_coordinates_byNucleotide_from3p_[CDS_i].append(partial_end_i)
			else:
				partial_end_i = self.CDS_end_[CDS_i] - n_partial_nucleotide_ + 1
				if partial_end_i <= self.CDS_start_[CDS_i]:
					self.CDS_partial_coordinates_byNucleotide_from3p_[CDS_i].append(self.CDS_start_[CDS_i])
				else:
					self.CDS_partial_coordinates_byNucleotide_from3p_[CDS_i].append(partial_end_i)
				self.CDS_partial_coordinates_byNucleotide_from3p_[CDS_i].append(self.CDS_end_[CDS_i])
			# print(CDS_i, n_partial_nucleotide_, self.CDS_start_[CDS_i], self.CDS_end_[CDS_i])
			# print(self.CDS_partial_coordinates_byNucleotide_from3p_[CDS_i])


def main():
	f_CDSannotation = open(sys.argv[1], 'r')
	partial_style = sys.argv[2]
	partial_length = sys.argv[3]

	CDS = LoadAnnotation(f_CDSannotation)
	CDS.get_CDSannotation()

	CDS_partial = GetPartialAnnotation(CDS.CDS_, CDS.CDS_start_, CDS.CDS_end_, CDS.CDS_strand_)
	if partial_style == 'portion':
		CDS_partial.get_partial_byPortion(int(partial_length))
		for k_CDS, v_coordinates in CDS_partial.CDS_partial_coordinates_byPortion_.items():
			print(k_CDS + "\t" + CDS.CDS_strand_[k_CDS] + "\t" + "\t".join([str(pos) for pos in v_coordinates]))
	if partial_style == 'nucleotide':
		
		####### splitting from either 5'end or 3'end 
		partial_direction = sys.argv[4]
		if partial_direction == '5p':
			CDS_partial.get_partial_byNucleotide_from5p(int(partial_length))
			for k_CDS, v_coordinates in CDS_partial.CDS_partial_coordinates_byNucleotide_from5p_.items():
				print(k_CDS + "\t" + "\t".join([str(pos) for pos in v_coordinates]) + "\t" + CDS.CDS_strand_[k_CDS])
		else:
			CDS_partial.get_partial_byNucleotide_from3p(int(partial_length))
			for k_CDS, v_coordinates in CDS_partial.CDS_partial_coordinates_byNucleotide_from3p_.items():
				print(k_CDS + "\t" + "\t".join([str(pos) for pos in v_coordinates]) + "\t" + CDS.CDS_strand_[k_CDS])


if __name__ == "__main__":
	main()
