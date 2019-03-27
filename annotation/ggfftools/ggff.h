#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "ggff.pb.h"
#include "leveldb/db.h"

namespace ggff {
void populate_subgraph(const std::string &ggff_col1, pb::SubGraph *sg);

pb::GGFF parse_ggff(const std::string &ggff_path);
void write_ggff(const pb::GGFF &g, std::ostream &oss = std::cout);
void write_subgraphs(const std::vector<pb::SubGraph> &graphs,
                     std::ostream &oss = std::cout);

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

  using HandleRecord = std::function<void(const pb::GGFFRecord &)>;
  void iterate(HandleRecord &handle_func,
               const std::string &record_key_prefix = "");

  std::vector<pb::SubGraph>
  intersect_of(const pb::SubGraph &g,
               const std::string &record_key_prefix = "");

  std::vector<pb::SubGraph>
  difference_of(const pb::SubGraph &g,
                const std::string &record_key_prefix = "");

  std::vector<pb::SubGraph> union_of(const pb::SubGraph &g,
                                     const std::string &record_key_prefix = "");

  std::vector<pb::SubGraph> xor_of(const pb::SubGraph &g,
                                   const std::string &record_key_prefix = "");

private:
  // {GFF3_TYPE}:{GFF3_SOURCE} -> GGFFRecord PROTO
  std::unique_ptr<leveldb::DB> index = nullptr;
};
} // namespace ggff
