// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Xiao Liang $
// $Authors: Xiao Liang, Alen Šarić $
// --------------------------------------------------------------------------
//

#pragma once

#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>

namespace OpenMS
{

  /**
      @ingroup Chemistry

       * @brief Constructs a DigestionEnzymeProtein from amino acid cleavage rules.
*
* @param name Name of the enzyme
* @param cut_before Set of amino acids before/after which cleavage occurs (e.g. "KR" for Trypsin)
* @param sense Whether cleavage is C-terminal or N-terminal
* @param nocut_after Set of amino acids that inhibit cleavage (e.g. "P" for Trypsin)
* @param synonyms Optional synonyms for the enzyme
* @param regex_description Optional description of the regex
*
* @throw Exception::MissingInformation if cut_before is empty
* @throw Exception::InvalidParameter if cut_before or nocut_after contain non-uppercase amino acid characters
*
* @note 'X' is automatically appended to cut_before to match any amino acid
  */
  class OPENMS_DLLAPI DigestionEnzymeProtein :
    public DigestionEnzyme
  {
  public:
    enum class Sense {C_TERM,N_TERM};
    /** @name Constructors
    */
    //@{

    /// Default constructor
    DigestionEnzymeProtein();

    /// Constructor from base class (adding defaults for the missing stuff)
    explicit DigestionEnzymeProtein(const DigestionEnzyme& d);

    /// Copy constructor
    DigestionEnzymeProtein(const DigestionEnzymeProtein&) = default;

    /// Move constructor
    DigestionEnzymeProtein(DigestionEnzymeProtein&&) = default;

    /// Detailed constructor
    explicit DigestionEnzymeProtein(const String& name,
                                    const String& cleavage_regex,
                                    const std::set<String>& synonyms = std::set<String>(),
                                    const String& regex_description = "",
                                    EmpiricalFormula n_term_gain = EmpiricalFormula("H"),
                                    EmpiricalFormula c_term_gain = EmpiricalFormula("OH"),
                                    const String& psi_id = "",
                                    const String& xtandem_id = "",
                                    Int comet_id = -1,
                                    Int msgf_id = -1,
                                    Int omssa_id = -1);

    explicit DigestionEnzymeProtein(const String& name,
                             const String& cut_before,
                             Sense sense,
                             const String& nocut_after = "",
                             const std::set<String>& synonyms = std::set<String>(),
                             const String& regex_description = "");

    ///  Destructor
    ~DigestionEnzymeProtein() override;
    //@}

    /** @name Assignment
     */
    //@{
    /// Assignment operator
    DigestionEnzymeProtein& operator=(const DigestionEnzymeProtein&) = default;

    /// Move assignment operator
    DigestionEnzymeProtein& operator=(DigestionEnzymeProtein&&) & = default;
    //@}

    /** Accessors
    */
    //@{
    /// sets the N-terminal gain
    void setNTermGain(const EmpiricalFormula& value);

    /// returns N-terminal gain
    EmpiricalFormula getNTermGain() const;

    /// sets the C-terminal gain
    void setCTermGain(const EmpiricalFormula& value);

    /// returns C-terminal gain
    EmpiricalFormula getCTermGain() const;

    /// sets the PSI ID
    void setPSIID(const String& value);

    /// returns the PSI ID
    String getPSIID() const;

    /// sets the X! Tandem enzyme ID
    void setXTandemID(const String& value);

    /// returns the X! Tandem enzyme ID
    String getXTandemID() const;

    /// sets the Comet enzyme ID
    void setCometID(Int value);

    /// returns the Comet enzyme ID
    Int getCometID() const;

    /// sets the MSGFPlus enzyme id
    void setMSGFID(Int value);

    /// returns the MSGFPlus enzyme id
    Int getMSGFID() const;

    /// sets the OMSSA enzyme ID
    void setOMSSAID(Int value);

    /// returns the OMSSA enzyme ID
    Int getOMSSAID() const;

    //@}

    /** @name Predicates
    */
    //@{
    /// equality operator
    bool operator==(const DigestionEnzymeProtein& enzyme) const;

    /// inequality operator
    bool operator!=(const DigestionEnzymeProtein& enzyme) const;

    // Note: comparison operator is not inherited. TODO rename and make virtual
    /// equality operator for regex
    bool operator==(const String& cleavage_regex) const;

    /// equality operator for regex
    bool operator!=(const String& cleavage_regex) const;

    /// order operator
    bool operator<(const DigestionEnzymeProtein& enzyme) const;
    //@}

    /**
       @brief Set the value of a member variable based on an entry from an input file

       Returns whether the key was recognized and the value set successfully.
    */
    bool setValueFromFile(const String& key, const String& value) override;

    /// ostream iterator to write the enzyme to a stream
    friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const DigestionEnzymeProtein& enzyme);

  protected:
    EmpiricalFormula n_term_gain_;

    EmpiricalFormula c_term_gain_;

    String psi_id_;

    String xtandem_id_;

    Int comet_id_;

    Int msgf_id_;

    Int omssa_id_;

    // @param cut_before: a set of Amino Acids, before which a cut in a given sequence should be set
    // @param nocut_after: a set of Amino Acids, which disvalidate a cut, even though a given Amino Acids from cut_before has been met
    // @param sense: the sense, as to how the sequence has to be read.
    String buildRegex_(String cut_before, const String& nocut_after,const DigestionEnzymeProtein::Sense& sense);
  };


  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const DigestionEnzymeProtein& enzyme);

  typedef DigestionEnzymeProtein Protease;
}
