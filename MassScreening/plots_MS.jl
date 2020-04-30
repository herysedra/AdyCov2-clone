using Plots,MAT, StatsPlots, Statistics,JLD2,CSV,DataFrames,Printf,ColorBrewer,Colors#, Images, ImageView,ImageDraw
    Plots.default(grid=false,legendfontsize=6, legendfontcolor=:gray30,titlefontsize=9, titlefontcolor=:gray30,xtickfontsize=9, ytickfontsize=9,tickfontcolor=:black,titlefontfamily="Cambria",xguidefontsize=9)
    StatsPlots.default(grid=false,legendfontsize=6, legendfontcolor=:gray30,titlefontsize=9, titlefontcolor=:gray30,xtickfontsize=9, ytickfontsize=9,tickfontcolor=:black,linewidth=false,xguidefontsize=9)
    colors=[:orange,:purple3,:maroon,:gold,:orangered,:grey,:purple,:ivory3,:chocolate1,:tan1,:rosybrown,:rosybrown2,:brown2,:brown3,:brown4,:deeppink3,:deeppink4]
        #colors2=ColorBrewer.palette("Pastel2",8);colors2=repeat(colors2,outer=[10])
        colors2=[ColorBrewer.palette("Set1",7);ColorBrewer.palette("Set2",7);ColorBrewer.palette("Set3",12);ColorBrewer.palette("Pastel2",8)]
        colors2_subset=[colors2[15],colors2[#=16=#26],colors2[9],colors2[#=15/25=#30],colors2[29],colors2[18]]
        colors2_subset2=[colors2[8],colors2[9],colors2[26],colors2[30]]#,colors2[12]]#,colors2[17],colors2[19]]
        colors2_subset3=ColorBrewer.palette("Set3",12);colors2_subset3=[colors2[1];colors2_subset3[4:8];colors2_subset3[10:12]]
        colors2_subset4=ColorBrewer.palette("Spectral",11)
        colors2_subset5=ColorBrewer.palette("RdYlGn",11);colors2_subset5=[colors2_subset5[1:5];reverse(colors2_subset5[6:10])]
        colors2_subset6=reverse(ColorBrewer.palette("YlOrRd",9));colors2_subset6[8]=ColorBrewer.palette("Set3",12)[1];colors2_subset6[9]=ColorBrewer.palette("Set3",12)[end]
        colors2_subset7=[ColorBrewer.palette("Set1",5)[1];ColorBrewer.palette("Set1",5)[5];reverse(ColorBrewer.palette("Accent",8))]
        c_R01_I=colors2_subset[2];c_R02_I=colors2_subset[3]
    markershapes=[:cross,:star4, :xcross, :star6, :hline]
    riskregionnames = ["Machakos/Muranga" "Mandera" "Baringo/Nakuru" "Nairobi" "Uasin Gishu" "Marsabit" "Garissa" "Nakuru/Narok" "Turkana" "Taita Taveta" "Kitui/Meru" "Kilifi/Mombasa" "Kericho/Kisumu" "Kilifi/Lamu" "Kakamega/Kisumu" "Wajir" "Kajiado/Kisumu" "Homa bay/Migori" "Samburu/Laikipia" "Kilifi/Kwale" "Total"]
    counties_names=["Baringo","Bomet","Bungoma","Busia","Elgeyo-Marakwet","Embu","Garissa","Homa Bay","Isiolo","Kajiado","Kakamega","Kericho","Kiambu","Kilifi","Kirinyaga","Kisii","Kisumu","Kitui","Kwale","Laikipia","Lamu","Machakos","Makueni","Mandera","Marsabit","Meru","Migori","Mombasa","Murang'a","Nairobi","Nakuru","Nandi","Narok","Nyamira","Nyandarua","Nyeri","Samburu","Siaya","Taita Taveta","Tana River","Tharaka-Nithi","Trans Nzoia","Turkana","Uasin Gishu","Vihiga","Wajir","West Pokot"]
        counties_names=[s[1:min(4,length(s))]   for s in counties_names]
        counties_names=["Barin","Bom","Bung","Busi","Elge","Emb","Garis","Hom","Isiol","Kajia","Kaka","Keric","Kiam","Kilifi","Kirin","Kisii","Kisu","Kitui","Kwal","Laiki","Lam","Mach","Maku","Mand","Mars","Meru","Migo","Mom","Mura","Nairo","Naku","Nand","Naro","Nya","Nyan","Nyer","Sam","Siay","Taita","Tana","Thar","Tran","Turk","Uasi","Vihig","Wajir","West"]
    wa_coords=[[300,450], [515,85], [165,360], [235,465], [115,300], [300,140], [510,380], [180, 430], [120, 130], [355, 615], [340, 375], [465, 630], [100, 380], [495, 530], [40, 355], [490, 255], [155, 495], [30, 440], [250, 290], [400, 670]]
    ages=[string(i)*"-"*string(i+4) for i=0:5:70];push!(ages,"75+")#"0-4" "5-9" "10-14" "15-19" "20-24" "25-29" "30-34" "35-39" "40-44" "45-49" ]

    S0=[4.138758e6, 867417.0, 2.326182e6, 8.084069e6, 3.229145e6, 459761.0, 999280.0, 1.979082e6, 926952.0, 340661.0, 2.381706e6, 2.126254e6, 2.960717e6, 786461.0, 7.478259e6, 781212.0, 2.114588e6, 4.094022e6, 569586.0, 917976.0]
    S012A=[279992, 264585, 246992, 206735, 220761, 211260, 181931, 130091, 111761, 82036, 55221, 42469, 34381, 23781, 16214, 18044]
    #@load "data/data_for_age_structuredmodel_16n_a.jld2" N_region_age;    S0perAge=[sum(N_region_age[:,a]) for a=1:16]; print(S0perAge) ;
    S0perAge=[5993115, 6202467, 6345902, 5285706, 4447468, 3854402, 3570565, 2650024, 2259166, 1786205, 1308564, 1118068, 869994, 658162, 514521, 697759]
    width1=400;length1=400;     width2=1000;length2=800;        ylims_boxplot=(6e5,Inf)
    Rotation=40;    xRotationV=90;   limit=1e2;

########### Verify scenarios
function curves1(folders,plot_folder)
    #curve_names=["cumIs","cumDeaths","I","Q","S","R"]
    curve_names=["cumIs"]
    if !isdir(plot_folder)   mkdir(plot_folder)   end
    for folder in folders
        data_files=readdir(folder)
        data_files=[s for s in data_files if endswith(s,".jld2")&&s!="stats.jld2"]
        for data_file in data_files
            plot_folder2=plot_folder*data_file[1:end-4]*"/"
            if !isdir(plot_folder2)   mkdir(plot_folder2)   end
            @load folder*data_file  sims_vector
            for wa=1:size(riskregionnames,2)-1
                p_curves=[plot(size=(width1,length1)) for i=1:1]#for i=1:6]
                for sim=1:size(sims_vector,1),i=1:1#6
                    d=[sims_vector[sim][i][t][wa][1] for t=1:366]
                    plot!(p_curves[i],d,color=colors2[i])
                end
                for i=1:1#=6=#   savefig(p_curves[i],plot_folder2*string(i)*"_"*curve_names[i]*"_wa"*string(wa)*"_"*replace(riskregionnames[wa],"/"=>"-")*".png")   end
            end
        end
    end
end

########## Mass Screening
function mysort2(folders)
    data_files=[]
    for folder in folders
        data_files2=readdir(folder);
        data_files2=[folder*s for s in data_files2 if endswith(s,".jld2")]
        data_files=[data_files ;data_files2]
    end
    return data_files
end
function make_plots_massScreening(folders,plot_folder)
    if !isdir(plot_folder)   mkdir(plot_folder)   end

    data_files=mysort2(folders)
    cumIs=[];  cumIs_4=[];
    xticks_labels=["1e4","2.5e4","5e4","7.5e4","1e5","2.5e5","5e5","7.5e5","1e6"]#repeat(["1e4","2.5e4","5e4","7.5e4","1e5","2.5e5","5e5","7.5e5","1e6"],outer=3);
    for i=1:size(data_files,1)
        @load data_files[i]  sims_vector;
        final_cum=[];   final_cum_4=[];
        for sim=1:size(sims_vector,1) #for each sim
            push!(final_cum,sum(sims_vector[sim][3][:,:,1:2])) # final cumulative I, summed for all ages
            push!(final_cum_4,sum(sims_vector[sim][3][4,:,1:2]))
        end
        push!(cumIs,median(final_cum));
        push!(cumIs_4,median(final_cum_4));

    end
    ##Plotting
    b1=bar(cumIs,title="Cases all Kenya",xticks=(1:10:size(cumIs,1),xticks_labels),color=colors2[1:10]);#display(b1)
        b2=bar(cumIs_4,title="Cases in Nairobi",xticks=(1:10:size(cumIs,1),xticks_labels),color=colors2[1:10]);#display(b2)


        #=cumIs_4_m5tom9=reshape([cumIs_4[(i-1)*10+j]   for i=1:9,j=5:9],9*5)
        b2_m5tom9=bar(cumIs_4_m5tom9,title="Cases in Nairobi",xticks=(1:5:size(cumIs,1),xticks_labels),color=colors2_subset6[5:9]);#display(b2_m5tom9)
        #scatter([cumIs_4[(i-1)*10+j]   for i=1:9,j=5:9],title="Cases in Nairobi",xticks=(1:5:size(cumIs,1),xticks_labels),color=repeat(colors2_subset6[5:9],outer=9))
        cumIs_4_m5tom9=[cumIs_4[(i-1)*10+j]   for i=1:9,j=5:9]
        plot(cumIs_4_m5tom9,title="Cases in Nairobi",xticks=(1:5:size(cumIs,1),xticks_labels),color=colors2_subset6[5:9],markershape=markershapes)
        #CSV.write(plot_folder*"statsCumIs_m5tom9.csv",DataFrame(cumIs_4_m5tom9))=#


        bb=bar(title="Cases in Kenya");n_sc=Int(floor(size(cumIs_4,1)/10))
        bb_4=bar(title="Cases in Nairobi");n_sc=Int(floor(size(cumIs_4,1)/10))
        #labels=["1e4","2.5e4","5e4","7.5e4","1e5","2.5e5","5e5","7.5e5","1e6"]
        colors2_subset7=ColorBrewer.palette("Set3",10);
        for i=1:n_sc
            bar!(bb_4,cumIs_4[(i-1)*10+1:(i-1)*10+10],color=colors2_subset6[i],label=xticks_labels[i])
            bar!(bb,cumIs[(i-1)*10+1:(i-1)*10+10],color=colors2_subset6[i],label=xticks_labels[i])
        end
        bar!(bb_4,xticks=(1:1:10,["m"*string(j) for j=1:10]),legend=:bottomright,xlabel="Mass screening performed during month X");#display(bb_4)
        bar!(bb,xticks=(1:1:10,["m"*string(j) for j=1:10]),legend=:bottomright,ylims=(5e6,Inf),xlabel="Mass screening performed during month X");#display(bb)

        #CSV.write(plot_folder*"statsCumIs.csv",DataFrame(x=repeat(xticks_labels,inner=10) .* repeat(["m"*string(j) for j=1:10],outer=9),cumIs=cumIs,cumIs_4=cumIs_4))

        savefig(b1,plot_folder*"massScreening_sc1_Kenya_1.png");savefig(b2,plot_folder*"massScreening_sc1_Nairobi_1.png");
        savefig(bb,plot_folder*"massScreening_sc1_Kenya_2.png");savefig(bb_4,plot_folder*"massScreening_sc1_Nairobi_2.png");
    return cumIs,cumIs_4
end
