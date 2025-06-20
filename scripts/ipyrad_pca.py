import ipyrad.analysis as ipa
import pandas as pd
import toyplot

data = "filtered_final_83.recode.snps.hdf5"
imap = {
    'lineataE': ['Jlin_03751LM', 'Jlin_03752LM', 'Jlin_03753LM', 'Jlin_03760RC', 'Jlin_03761RC', 'Jlin_03762RC', 'Jlin_08228Rc', 'Jlin_08230Rc', 'Jlin_08231Rc', 'Jmul_06621PT', 'Jmul_06623PT', 'Jmul_06625PT', 'Jmul_08311LD', 'Jmul_08314LD', 'Jmul_08315LD', 'Jmul_08316Ld', 'Jmul_08317Ld', 'Jmul_08318Ld'], 
    'luxata': ['Jlux_08945BO', 'Jlux_08946BO', 'Jlux_08947BO', 'Jlux_17668BF', 'Jlux_17669BF', 'Jlux_17670BF'], 
    'lineataW': ['Jmul_03771SR', 'Jmul_03772SR', 'Jmul_03773SR', 'Jmul_03778CE', 'Jmul_03780CE', 'Jmul_03781CE', 'Jmul_03795AR', 'Jmul_03796ar', 'Jmul_03797AR', 'Jmul_03851CE', 'Jmul_03852CE', 'Jmul_03853CE', 'Jmul_06589SO', 'Jmul_06594SO', 'Jmul_06604SO', 'Jmul_08352SO', 'Jmul_08353CE', 'Jmul_08354SO', 'Jmul_08377SO', 'Jmul_08400SO', 'Jspp_06617SU', 'Jspp_06619SU', 'Jspp_06620SU'], 
    'darwini': ['Jmul_08671LP', 'Jmul_08672LP', 'Jmul_08674LP', 'Jmul_08721LC', 'Jmul_08722LC', 'Jmul_08723LC', 'Jmul_08748LL', 'Jmul_08749LL', 'Jmul_08750LL', 'Jmul_09891LP', 'Jmul_09892LP', 'Jmul_09893LP'], 
    'lineataN': ['Jmul_17901MT', 'Jmul_17902MT', 'Jmul_17903MT', 'Jmul_17909RS', 'Jmul_17912RS', 'Jmul_17915RS', 'Jspp_06607SU', 'Jspp_06616SU', 'Jspp_06618SU'], 
    'onca': ['Jonc_03621RN', 'Jonc_03626OG', 'Jonc_03627RN', 'Jonc_03696RY', 'Jonc_03697RY', 'Jonc_03699RY', 'Jonc_03716T1', 'Jonc_03717T1', 'Jonc_03718T1', 'Jonc_03734T3', 'Jonc_03735T3', 'Jonc_03736T3', 'Jonc_03739T4', 'Jonc_03740T4', 'Jonc_03741T4']
}


pca = ipa.pca(
    data=data,
    imap=imap,
    impute_method="sample",
)


pca.run(nreplicates=100, seed=12345)

# store the PC axes as a dataframe
df = pd.DataFrame(pca.pcaxes[0], index=pca.names)

# write the PC axes to a CSV file
df.to_csv("pca.csv")

pca.draw(0, 2);

pca.draw(outfile="mypca.pdf")
