#pragma once

#include <memory>

#include "dgt_defines.hpp"

namespace dgt {

class When
{
  public:
    When() = default;
    virtual ~When();
    virtual void limit_time(real&) const {}
    virtual bool now(int const step, real const time) = 0;
};

class AtAlways : public When
{
  public:
    bool now(int const step, real const time) override;
};

class AtNever : public When
{
  public:
    bool now(int const step, real const time) override;
};

class AtStep : public When
{
  protected:
    int m_step;
  public:
    AtStep(int const step);
    bool now(int const step, real const time) override;
};

class AtStepPeriodically : public When
{
  protected:
    int m_frequency;
    int m_last_step;
  public:
    AtStepPeriodically(int const frequency);
    bool now(int const step, real const time) override;
};

class AtTime : public When
{
  protected:
    real m_time;
    bool m_has_occurred;
  public:
    AtTime(real const time);
    bool now(int const step, real const time) override;
};

class AtExactTime : public AtTime
{
  public:
    AtExactTime(real const time);
    void limit_time(real& new_time) const override;
};

class AtTimePeriodically : public When
{
  protected:
    real m_frequency;
    int m_iteration;
  public:
    AtTimePeriodically(real const frequency);
    bool now(int const step, real const time) override;
};

class AtExactTimePeriodically : public AtTimePeriodically
{
  public:
    AtExactTimePeriodically(real const frequency);
    void limit_time(real& new_time) const override;
};

using WhenPtr = std::shared_ptr<When>;

class AtEither : public When
{
  protected:
    WhenPtr m_first = nullptr;
    WhenPtr m_second = nullptr;
  public:
    AtEither(WhenPtr first, WhenPtr second);
    void limit_time(real& new_time) const override;
    bool now(int const step, real const time) override;
};

WhenPtr combine_either(WhenPtr a, WhenPtr b);

}
