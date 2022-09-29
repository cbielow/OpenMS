// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: $
// $Authors: Moritz Berger, Tetana Krymovska$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CONCEPT/Colorizer.h>
///////////////////////////

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

using namespace OpenMS;
using namespace std;


template <typename T>
//convert any input to string
string convertToString ( T var_input )
{
    ostringstream ss;
    ss << var_input;
    return ss.str();
}

START_TEST(Colorizer(),"$Id$")

 //Test variables
char test_char = 'a';
int test_int = 15;
float test_float = 2094.5892;
string test_string = " !#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

//ANSI codes

#ifdef OPENMS_WINDOWSPLATFORM

    string const blackANSI        = "";
    string const redANSI          = "";
    string const greenANSI        = "";
    string const yellowANSI       = "";
    string const blueANSI         = "";
    string const magentaANSI      = "";
    string const cyanANSI         = "";
    string const whiteANSI        = "";
    string const resetColorANSI   = "";

#elif defined(__linux__) || defined(__OSX__)

    string const blackANSI        = "\e[30m";
    string const redANSI          = "\e[31m";
    string const greenANSI        = "\e[32m";
    string const yellowANSI       = "\e[33m";
    string const blueANSI         = "\e[34m";
    string const magentaANSI      = "\e[35m";
    string const cyanANSI         = "\e[36m";;
    string const whiteANSI        = "\e[37m";
    string const resetColorANSI   = "\e[0m";

#endif


START_SECTION(Colorizer::colorStream(ostream& stream) const)
{
    //without text
    stringstream test_stream;
    ColorizerTester c(Color::BLACK);

    c.colorStreamSimple(test_stream);
    TEST_EQUAL(test_stream.str(), blackANSI)

    //with text
    test_stream.str(string());
    test_stream.clear();

    test_stream << c(test_string);
    c.colorStreamSimple(test_stream);
    TEST_EQUAL(test_stream.str(),test_string+blackANSI)
}
END_SECTION


START_SECTION(Colorizer::outputToStream(ostream& o_stream))
{
    //without text
    stringstream test_stream;
    ColorizerTester c(Color::CYAN);
    //c();
    c.outputToStreamSimple(test_stream);
    TEST_EQUAL(test_stream.str(), cyanANSI+resetColorANSI)
    

    // with text
    test_stream.str(string());
    test_stream.clear();

    test_stream << c(test_string);
    c.outputToStreamSimple(test_stream);

    TEST_EQUAL(test_stream.str(),test_string+cyanANSI+test_string+resetColorANSI) //modified TEST_EQUAL
}
END_SECTION

START_SECTION(Colorizer::resetColor(ostream& stream))
{
    stringstream test_stream;
    ColorizerTester c(Color::GREEN);

    test_stream << c(test_string);
    c.resetColorSimple(test_stream);
    TEST_EQUAL(test_stream.str(), test_string+resetColorANSI)
}
END_SECTION

START_SECTION(Colorizer::getDataAsString())
{
    stringstream test_stream;
    ColorizerTester c(Color::RED);

    test_stream << c(test_string);
    test_stream << c.getDataAsString_();
    TEST_EQUAL(test_stream.str(), test_string+test_string)
}
END_SECTION

// START_SECTION(Colorizer::reset())
// {
//     stringstream test_stream;
//     test_stream << green() 
//         << "green text" << 89 << "$" << " " 
//         << green.reset() << "default text" << red() << 11 << red.reset() << "A";
//     TEST_EQUAL(test_stream.str(),greenANSI
//                                 +"green text89$ "
//                                 +greenANSI+resetColorANSI
//                                 +"default text"
//                                 +redANSI+"11"+redANSI+resetColorANSI+"A")
// }
// END_SECTION

START_SECTION("Testing Colorizer instances")
{
    //Check that the colorized input contains the original text and according ASCI codes

    stringstream colored_stream;

    colored_stream << black(test_string);
    TEST_EQUAL(colored_stream.str(), test_string)
    colored_stream.str(string());
    colored_stream.clear();

    // TEST_EQUAL(colored_stream.str().find(blackANSI) != std::string::npos,1) //no
    // TEST_EQUAL(colored_stream.str().find(resetColorANSI) != std::string::npos,1) //no
    // TEST_EQUAL(colored_stream.str().find(test_string) != std::string::npos,1) //yes delete

    colored_stream << red(test_string);
    TEST_EQUAL(colored_stream.str(), test_string)
    colored_stream.str(string());
    colored_stream.clear();

    colored_stream << green(test_string);
    TEST_EQUAL(colored_stream.str(), test_string)
    colored_stream.str(string());
    colored_stream.clear();

    colored_stream << yellow(test_string);
    TEST_EQUAL(colored_stream.str(), test_string)
    colored_stream.str(string());
    colored_stream.clear();

    colored_stream << blue(test_string);
    TEST_EQUAL(colored_stream.str(), test_string)
    colored_stream.str(string());
    colored_stream.clear();

    colored_stream << magenta(test_string);
    TEST_EQUAL(colored_stream.str(), test_string)
    colored_stream.str(string());
    colored_stream.clear();

    colored_stream << cyan(test_string);
    TEST_EQUAL(colored_stream.str(), test_string)
    colored_stream.str(string());
    colored_stream.clear();

    colored_stream << white(test_string);
    TEST_EQUAL(colored_stream.str(), test_string)
    colored_stream.str(string());
    colored_stream.clear();
}
END_SECTION

//testing various inputs for colorizing
START_SECTION("Testing Colorizer inputs")
{
    stringstream colored_stream;
    string comparison_string;

    //char////////////////////////////////////
    colored_stream << yellow(test_char);
    comparison_string.append(convertToString(test_char));
    TEST_EQUAL(colored_stream.str(), comparison_string)

    //clearing streams
    colored_stream.str(string());
    colored_stream.clear();
    comparison_string = "";

    //int///////////////////////////////
    colored_stream << green(test_int);
    comparison_string.append(convertToString(test_int));
    TEST_EQUAL(colored_stream.str(), comparison_string)

    //clearing streams
    colored_stream.str(string());
    colored_stream.clear();
    comparison_string = "";

    //float///////////////////////////////
    colored_stream << red(test_float);
    comparison_string.append(convertToString(test_float));
    TEST_EQUAL(colored_stream.str(), comparison_string)

    //clearing streams
    colored_stream.str(string());
    colored_stream.clear();
    comparison_string = "";
}
END_SECTION

START_SECTION(Colorizer& operator()())
{
    stringstream test_stream;
    test_stream << green() 
                << "green text" 
                << 123 << "!" 
                << " ";
    TEST_EQUAL(test_stream.str(),"green text123! ")
}
END_SECTION

END_TEST