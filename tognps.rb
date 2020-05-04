# translates a project space into a GNPS input

Compound = Struct.new(:name, :mz, :rt, :quant)

compounds = []
#>quantification "SPF3_Thymus_GH8_01_22419":(7532.0);"SPF4_Jeju5_GD2_01_22583":(6756.0);"GF4_Jeju4_BD1_01_22469":(7024.0);"SPF4_Kid3_GF1_01_22610":(7220.0);"Blank_GA1_01_22653":(9120.0);"SPF4_Heart4_GG6_01_22629":(7004.0);"SPF4_Spleen3_GB2_01_22555":(6776.0);"Blank_GA1_01_22510":(8076.0);"GF4_Brain2_BH10_01_22534":(6684.0);"SPF3_Lung3_GG9_01_22406":(6692.0);"Blank_GA1_01_22609":(9748.0)



File.open("gnps.mgf","w") {|io|
  Dir.glob("#{ARGV[0]}/*/spectrum.ms").sort_by {|x| File.basename(File.dirname(x)).split("_").last.to_i}.each {|x|
    spec = File.readlines(x)
    name = spec.grep(/>compound/).first.chomp.split(/\s+/,2).last.strip
    mz = spec.grep(/>parentmass (.+)/).first.chomp.split(/\s+/,2).last.to_f
    rt = spec.grep(/>rt (.+)/).first.chomp.split(/\s+/,2).last[0...-1].to_f
    quant = Hash[spec.grep(/>quantification/).first.chomp.split(/\s+/,2).last.split(";").map {|x| x =~ /"([^"]+)":\(([^\)]+)\)/; [$1,Float($2)]}]
    c = Compound.new(name,mz,rt,quant)
    io.puts "BEGIN IONS"
    io.puts "FEATURE_ID=#{c.name}"
    io.puts "PEPMASS=#{mz}"
    io.puts "CHARGE=1+"
    io.puts "MSLEVEL=2"
    io.puts "RTINSECONDS=#{rt}"
    io.puts "SCANS=#{name}"
    spec.slice_before(/^>(ms1|ms2|coll)/).drop(1).each {|batch|
      if batch[0]=~/^>(ms2|col)/
        batch[1..-1].grep(/^\d/).each {|peak|
          io.puts("#{peak}")
        }
      end
    }
    io.puts("END IONS\n")
    compounds<<c
  } 
};nil

# write quant
require "csv"
samples = compounds.flat_map {|x| x.quant.keys}.uniq.sort;nil
CSV.open("gnps_quant.csv","w") {|io|
  #row ID,row m/z,row retention time,correlation group ID,annotation network number,best ion,auto MS2 verify,identified by n=,partners,neutral M mass,10765_P4_RE9_01_482.mzXML Peak area,11035_P4_RB4_01_431.mzXML Peak area,10715_P4_RA4_01_415.mzXML Peak area,10712_P4_RH3_01_521.mzXML Peak area,14153_P4_RH8_01_526.mzXML Peak area,18404_P3_RF6_01_370.mzXML Peak area,15826_P2_RF6_01_238.mzXML Peak area,14880_P3_RA3_01_290.mzXML Peak area,13917_P4_RA10_01_423.mzXML Peak area,32754_P5_RE6_01_607.mzXML Peak area,1820_P3_RB5_01_309.mzXML Peak area,29342_P5_RF3_01_620.mzXML Peak area,33428_P6_RB2_01_678.mzXML Peak area,32764_P3_RF12_01_377.mzXML Peak area,32872_P2_RH11_01_276.mzXML
  row = ['row ID', 'row m/z', 'row retention time']
  samples.each {|s| row << s}

  io << row

  compounds.each {|c|
    row = [c.name, c.mz, c.rt]
    samples.each {|s|
      row << (c.quant[s] || 0.0)
    }

    io << row
  }



};nil
