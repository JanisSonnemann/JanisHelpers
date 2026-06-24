# FlowJo WSP XML structure — verified against tests/fixtures/minimal.wsp
#
# Root sample tree node name (child of <Sample>):          SampleNode
# Population container element:                            Subpopulations
# Individual population element:                           Population
# Population count attribute:                              count
# Population name attribute:                               name
# Keyword container element (child of <Sample>):           Keywords (capital K)
# Keyword child element:                                   Keyword
# Keyword name/value attributes:                           name, value
# Statistic element name:                                  Statistic
# Statistic type attribute (Median/Mean/etc.):             name
# Statistic channel attribute:                             id  (NOT "channel")
# Statistic value attribute:                               value
# NOTE: Statistic nodes are siblings of Population inside Subpopulations, NOT children of Population
# SampleRef ID attribute:                                  sampleID
# DataSet sample ID attribute:                             sampleID
# Namespace prefix required for XPath:                     none (xml_ns_strip not needed)
