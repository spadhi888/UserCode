/** \class DoubleValueProducer 
 *
 * \author Fedor Ratnikov (inpired by FWCore/Framework/test/stubs/ToyModules.cc)
 *
 * \version $Id: DoubleValueProducer.cc,v 1.1 2009/02/09 17:45:47 fedor Exp $
 *
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

class DoubleValueProducer : public edm::EDProducer {
 public:
    explicit DoubleValueProducer(edm::ParameterSet const& p) : 
      mValue (p.getUntrackedParameter<double>("value")) {
      produces<double>();
    }
    virtual ~DoubleValueProducer() { }
  virtual void produce(edm::Event& e, edm::EventSetup const& c) {
    std::auto_ptr<double> p(new double (mValue));
    e.put(p);
  }
 private:
  double mValue;
};

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( DoubleValueProducer );


