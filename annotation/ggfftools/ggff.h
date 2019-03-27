#include <memory>
#include <string>

#include "ggff.pb.h"
#include "leveldb/db.h"

namespace ggff {
void populate_subgraph(const std::string &ggff_col1, pb::SubGraph *sg);

pb::GGFF parse_ggff(const std::string &ggff_path);

std::string record_key(const pb::GGFFRecord &rec);

pb::Strand operator&(const pb::Strand &s1, const pb::Strand &s2);
pb::Strand operator-(const pb::Strand &s1, const pb::Strand &s2);
pb::Strand operator|(const pb::Strand &s1, const pb::Strand &s2);
pb::Strand operator^(const pb::Strand &s1, const pb::Strand &s2);

pb::NodeSlice operator&(const pb::NodeSlice &sl1, const pb::NodeSlice &sl2);
pb::NodeSlice operator-(const pb::NodeSlice &sl1, const pb::NodeSlice &sl2);
pb::NodeSlice operator|(const pb::NodeSlice &sl1, const pb::NodeSlice &sl2);
pb::NodeSlice operator^(const pb::NodeSlice &sl1, const pb::NodeSlice &sl2);

pb::SubGraph operator&(const pb::SubGraph &sg1, const pb::SubGraph &sg2);
pb::SubGraph operator-(const pb::SubGraph &sg1, const pb::SubGraph &sg2);
pb::SubGraph operator|(const pb::SubGraph &sg1, const pb::SubGraph &sg2);
pb::SubGraph operator^(const pb::SubGraph &sg1, const pb::SubGraph &sg2);

class GGFFIndex {
public:
  explicit GGFFIndex(const std::string &ggff_index_path);
  GGFFIndex() = delete;

  void add_ggff(const std::string &ggff_path);

private:
  // {GFF3_TYPE}:{GFF3_SOURCE} -> GGFFRecord PROTO
  std::unique_ptr<leveldb::DB> index = nullptr;
};
} // namespace ggff
