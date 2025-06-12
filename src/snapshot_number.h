#ifndef SNAPSHOT_NUMBER_H_INCLUDED
#define SNAPSHOT_NUMBER_H_INCLUDED

#include <assert.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include "config_parser.h"
#include "datatypes.h"

class SnapshotNumber_t
{
protected:
  int SnapshotIndex;
  int SnapshotId;

public:
  SnapshotNumber_t()
  {
    SnapshotIndex = SpecialConst::NullSnapshotId;
    SnapshotId = SpecialConst::NullSnapshotId;
  }
  SnapshotNumber_t(SnapshotNumber_t &sn) : SnapshotId(sn.SnapshotId), SnapshotIndex(sn.SnapshotIndex)
  {
  }
  SnapshotNumber_t &operator=(SnapshotNumber_t &sn)
  {
    SnapshotIndex = sn.SnapshotIndex;
    SnapshotId = sn.SnapshotId;
    return *this;
  }
  void ResetSnapshotNumber()
  { // reset is not destructon! when destructor is called, the data content no matter matters.
    SnapshotIndex = SpecialConst::NullSnapshotId;
    SnapshotId = SpecialConst::NullSnapshotId;
  }
  void FormatSnapshotId(std::stringstream &ss);
  void SetSnapshotIndex(int snapshot_index);
  int GetSnapshotIndex() const;
  int GetSnapshotId() const;
  int GetSnapshotId(const int &snapshot_index) const;
};

inline int SnapshotNumber_t::GetSnapshotIndex() const
{
  return SnapshotIndex;
}

inline int SnapshotNumber_t::GetSnapshotId() const
{
  return SnapshotId;
}

/* Used to set SnapshotId member, and to retrieve the snapshot number of a given
 * snapshot index in different parts of the code. */
inline int SnapshotNumber_t::GetSnapshotId(const int &snapshot_index) const
{
  if (HBTConfig.SnapshotIdList.empty())
    return snapshot_index;
  else
    return HBTConfig.SnapshotIdList[snapshot_index];
}

inline void SnapshotNumber_t::FormatSnapshotId(stringstream &ss)
{
  ss << std::setw(3) << std::setfill('0') << SnapshotId;
}

inline void SnapshotNumber_t::SetSnapshotIndex(int snapshot_index)
{
  assert(snapshot_index >= HBTConfig.MinSnapshotIndex && snapshot_index <= HBTConfig.MaxSnapshotIndex);

  SnapshotIndex = snapshot_index;
  SnapshotId = GetSnapshotId(snapshot_index);
}

#endif