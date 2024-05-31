
import sys

class LoadAnnotation:
	"""load CDS annotation and TSS category"""

	def __init__(self, filename_CDSannotation, filename_tssCategory):
		self.filename_CDSannotation_ = filename_CDSannotation
		self.filename_tssCategory_ = filename_tssCategory
		self.CDS_ = []
		self.CDS_start_ = {}
		self.CDS_end_ = {}
		self.CDS_strand_ = {}
		self.iTSS_coordinate_ = []
		self.iTSS_strand_ = {}

	def get_CDSannotation(self):
		for line in self.filename_CDSannotation_:
			line = line.rstrip().split("\t")
			CDS = line[0]
			self.CDS_.append(CDS)
			self.CDS_start_[CDS] = int(line[1])
			self.CDS_end_[CDS] = int(line[2])
			self.CDS_strand_[CDS] = line[3]

	def get_tssCategory(self):
		next(self.filename_tssCategory_)
		for line in self.filename_tssCategory_:
			line = line.rstrip().split(",")
			tss_coordinate = int(line[0])
			if line[3] == 'YES':
				self.iTSS_coordinate_.append(tss_coordinate)
				self.iTSS_strand_[tss_coordinate] = line[1]


class GetTSSinCDS:
	"""assign iTSS to CDSs. get number of iTSS in each CDS"""

	def __init__(self, list_CDS, dict_CDS_start, dict_CDS_end, dict_CDS_strand, list_iTSS, dict_iTSS_strand):
		self.CDS_ = list_CDS
		self.CDS_start_ = dict_CDS_start
		self.CDS_end_ = dict_CDS_end
		self.CDS_strand_ = dict_CDS_strand
		self.iTSS_coordinate_ = list_iTSS
		self.iTSS_strand_ = dict_iTSS_strand
		self.CDS_iTSS_ = {}

	def get_iTSS_inCDS(self):
		for CDS_i in self.CDS_:
			for iTSS_j in self.iTSS_coordinate_:
				if self.iTSS_strand_[iTSS_j] == self.CDS_strand_[CDS_i]:
					if iTSS_j <= self.CDS_end_[CDS_i]:
						if iTSS_j >= self.CDS_start_[CDS_i]:
							# print(CDS_i, iTSS_j)
							self.CDS_iTSS_[CDS_i] = self.CDS_iTSS_.get(CDS_i, 0) + 1


def main():
	f_CDSannotation = open(sys.argv[1], 'r')
	f_tssCategory = open(sys.argv[2], 'r')

	CDS = LoadAnnotation(f_CDSannotation, f_tssCategory)
	CDS.get_CDSannotation()
	CDS.get_tssCategory()

	############## test session
	# print(CDS.iTSS_coordinate_)
	# print(len(CDS.iTSS_coordinate_))

	# print(CDS.iTSS_strand_)
	# print(len(CDS.iTSS_strand_))
	############## end of test session

	iTSS_count = GetTSSinCDS(CDS.CDS_, CDS.CDS_start_, CDS.CDS_end_, CDS.CDS_strand_, CDS.iTSS_coordinate_, CDS.iTSS_strand_)
	iTSS_count.get_iTSS_inCDS()

	############## test session
	# print(iTSS_count.CDS_iTSS_)
	# print(len(iTSS_count.CDS_iTSS_))
	############## end of test session

	for CDS_i, iTSS_j in iTSS_count.CDS_iTSS_.items():
		print(CDS_i + "\t" + str(iTSS_j))


if __name__ == "__main__":
	main()