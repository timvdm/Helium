#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/depict/depict.h"
#include "../../src/depict/svgpainter.h"

using Helium::Chemist::Molecule;
using namespace boost::python;

// strange bug on win32... :(
#ifdef WIN32
std::stringstream m_os;
#endif

class SVGPainter : public Helium::SVGPainter
{
  public:
    SVGPainter() : Helium::SVGPainter(m_os)
    {
      m_os.str("");
    }

    std::string output()
    {
      return m_os.str() + "</svg>";
    }

#ifndef WIN32
  private:
    std::stringstream m_os;
#endif
};

bool draw_molecule(Helium::Depict &self, const Molecule &mol, const Helium::RingSet<Molecule> &rings, const list &origCoords)
{
  std::vector<std::pair<double, double> > coords;
  for (std::size_t i = 0; i < len(origCoords); ++i) {
    tuple t = extract<tuple>(origCoords[i]);
    double x = extract<double>(t[0]);
    double y = extract<double>(t[1]);
    coords.push_back(std::make_pair(x, y));
  }
  return self.drawMolecule(mol, rings, coords);
}

Helium::Depict* Depict_ctor(SVGPainter *painter)
{
  return new Helium::Depict(painter);
}

void export_depict()
{

  class_<SVGPainter, boost::noncopyable>("SVGPainter", init<>())
    .def("output", &SVGPainter::output)
    ;

  scope in_Depict = class_<Helium::Depict, boost::noncopyable>("Depict", no_init)
    //.def(init<SVGPainter*>())
    .def("__init__", make_constructor(&Depict_ctor))
    .def("drawMolecule", &draw_molecule)
    .def("setBondLength", &Helium::Depict::setBondLength)
    .def("bondLength", &Helium::Depict::bondLength)
    .def("setPenWidth", &Helium::Depict::setPenWidth)
    .def("penWidth", &Helium::Depict::penWidth)
    .def("setBondSpacing", &Helium::Depict::setBondSpacing)
    .def("bondSpacing", &Helium::Depict::bondSpacing)
    .def("setBondWidth", &Helium::Depict::setBondWidth)
    .def("bondWidth", &Helium::Depict::bondWidth)
    .def("setOption", &Helium::Depict::setOption)
    .def("options", &Helium::Depict::options)
    .def("clearOptions", &Helium::Depict::clearOptions)
    .def("setFontFamily", &Helium::Depict::setFontFamily)
    .def("fontFamily", &Helium::Depict::fontFamily, return_value_policy<copy_const_reference>())
    .def("setFontSize", &Helium::Depict::setFontSize)
    .def("fontSize", &Helium::Depict::fontSize)
    .def("setBondColor", &Helium::Depict::setBondColor)
    ;

  enum_<Helium::Depict::Options>("Options")
    .value("BlackWhiteAtoms", Helium::Depict::BlackWhiteAtoms)
    .value("NoMargin", Helium::Depict::NoMargin)
    .value("DrawTermC", Helium::Depict::DrawTermC)
    .value("DrawAllC", Helium::Depict::DrawAllC)
    .value("NoWedgeHashGen", Helium::Depict::NoWedgeHashGen)
    .value("AsymmetricDoubleBond", Helium::Depict::AsymmetricDoubleBond)
    .value("AromaticCircle", Helium::Depict::AromaticCircle)
    .value("AromaticHash", Helium::Depict::AromaticHash)
    ;
}
