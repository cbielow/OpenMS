// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow, Moritz Berger, Tetana Krymovska$
// $Authors: Chris Bielow, Moritz Berger, Tetana Krymovska$
// --------------------------------------------------------------------------

#pragma once

//#include <OpenMS/KERNEL/MSExperiment.h>

#include <array>
#include <iosfwd>
#include <sstream>
#include <OpenMS/config.h>

//////
#include <unistd.h>
#include <stdio.h>
/////


namespace OpenMS
{
  /// Text colors for console output
  enum class Color
  { // internal note: the order is important here! See Colorizer::colors_
    BLACK,
    RED,
    GREEN,
    YELLOW,
    BLUE,
    MAGENTA,
    CYAN,
    WHITE,
    RESET, ///< reset the color to the previous (default?) color
  };

  /**
   * @brief A class, that provides options for colored output with the "<<" operator for output streams (cout, cerr)
   *
   */

  class OPENMS_DLLAPI Colorizer
  {
    
friend class ConsoleUtils;

public:
    /// Constructor
    Colorizer(const Color color);

    /// Constructor
    Colorizer();

    // Copy constructor
    Colorizer(const Colorizer &rhs);

    /// Destructor
    ~Colorizer();

    /// insetrion Operator
    friend std::ostream& operator<<(std::ostream& o_stream, Colorizer& col);


    /// Bracket Operator
    Colorizer& operator()()
    {
      reset_ = false;
      this->input_.str(""); // clear the stream
      return *this;
    }

    /// Bracket Operator
    template<typename T>
    Colorizer& operator()(T s)
    {
      this->input_.str(""); // clear the stringstream
      this->input_ << s;    // add new data
      reset_ = true;
      return *this;
    }

     ///
    Colorizer& reset();


protected:

    const int color_;
    ///
    void outputToStream(std::ostream& o_stream);

    ///
    void colorStream(std::ostream& stream) const;

    ///
    void resetColor(std::ostream& stream);

    ///
    std::string getDataAsString();

    bool getReset()
    {
      return this->reset_;
    }

    std::string getInput()
    {
        std::stringstream input_out;
        input_out << input_.rdbuf();

      return input_out.str();
     }

    const char* getColor_()
    {
      return this->colors_[color_];
    }

    const char* getResetColor_(){
      return this->colors_[int(Color::RESET)];
    }

private:
    // const int color_;

    /// input in Colorizer object to be colored
    std::stringstream input_;

    /// 
    bool reset_ = true;


/**
 * @brief constant string array which saves the Linux color codes.
 * 0=black
 * 1=red
 * 2=green
 * 3=yellow
 * 4=blue
 * 5=magenta
 * 6=cyan
 * 7=white
 * 8=default console color (reset)
 *
 */
#if defined(_WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(_WIN64)
    /// 
    inline static constexpr std::array<const int, 9> colors_ {16, 12, 10, 14, 9, 13, 11, 15, 15};

#elif defined(__linux__) || defined(__OSX__)
    /// 
    inline static constexpr std::array<const char*, 9> colors_ {"\033[30m", "\033[31m", "\033[32m", "\033[33m", "\033[34m", "\033[35m", "\033[36m", "\033[37m", "\033[0m"};

#endif
  };

  // declaration of all colorizer object.
  extern OPENMS_DLLAPI Colorizer black;
  extern OPENMS_DLLAPI Colorizer red;
  extern OPENMS_DLLAPI Colorizer green;
  extern OPENMS_DLLAPI Colorizer yellow;
  extern OPENMS_DLLAPI Colorizer blue;
  extern OPENMS_DLLAPI Colorizer magenta;
  extern OPENMS_DLLAPI Colorizer cyan;
  extern OPENMS_DLLAPI Colorizer white;
  // extern OPENMS_DLLAPI Colorizer reset_color;   ///< reset the color to default, alias for 'make_default_color'
  //extern /*OPENMS_DLLAPI*/ Colorizer default_color; ///< reset the color to default, alias for 'reset_color'
  
  //Stream operator declaration
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& o_stream, 
                                        OpenMS::Colorizer& col);




/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
class OPENMS_DLLAPI ColorizerTester: public Colorizer

      /*only used in Colorizer_test.cpp for testing Colorizer method
      functionality. While using Colorizer instances, only cout or cerr
      streams are being colorized. 

      Methods with underscore "_" at the name end are used to
      reveal functionality of protected class functions in the unit test.
      
      For separate Colorizer method testing, stringstreams are used instead. 
      This class modified these methods to use Colorizer functionality on 
      stringstreams instead of cout/cerr streams. Such methods have "Simple"
      attached to the method name.
      
      This class should not be used outside Colorizer_test.cpp
      */
{
    public:

    ///Constructor
    // ColorizerTester(const Color color);
    ColorizerTester(const Color color) : Colorizer(color){}; 


    /// Default destructor
    ~ColorizerTester(){};

    //original methods revealing functionality for public use
    void outputToStream_(std::ostream& o_stream){this->outputToStream(o_stream);}

    void colorStream_(std::ostream& stream) const{this->colorStream(stream);}

    void resetColor_(std::ostream& stream){this->resetColor(stream);}

    std::string getDataAsString_(){return this->getDataAsString();}


    //modified methods for stringstream instead of cerr/cout
    void outputToStreamSimple(std::ostream& o_stream);

    void colorStreamSimple(std::ostream& stream);

    void resetColorSimple(std::ostream& stream);

   // bool getReset_();

    // std::string getDataAsStringSimple(); //not needed
};

}
