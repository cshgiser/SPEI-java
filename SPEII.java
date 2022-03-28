package myPackage;
import java.lang.Math.*;
import org.apache.commons.math3.*;
import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.special.Gamma;

public class SPEII {
    private double[][] temp;  //年 月
    private double[][] pre;   //年 月
    private double[][] mothlyD;  //年 月
    private int k=1;
    private final int yearnum;
    private final int monthnum;
    private double latitude;
    private final double[][] xacc;

    public SPEII(double[][] temp,double[][] pre,double latitude,int scale,int startYear){
        this.k=scale;
        this.pre = pre;
        this.temp = temp;
        yearnum = temp.length;
        monthnum = temp[0].length;
        mothlyD = new double[yearnum][monthnum];  // monthly difference between precipitation and PET

        YearMonthlyDI ymDI;
        ymDI = new YearMonthlyDI(temp,pre,latitude,startYear);
        mothlyD=ymDI.getResult();

        this.xacc = calculateX();
    }

    public double[][] getResult(){
        return calculateSPEI();
    }

    private double[][] calculateSPEI(){
        double[] seasonSeries = new double[yearnum];
        double[][] spei = new double[yearnum][monthnum];

        for (int j = 0;j< monthnum;j++){
            // 周期 一年 12 months
            int kk=0;
            for(int i=0;i<yearnum;i++)
            {
                seasonSeries[kk] = xacc[i][j];
                kk++;
            }
            int nSeason = kk;
            double[] seasonXaccsort = bubbleSort(seasonSeries);
            double w0=0;
            double w1=0;
            double w2=0;
            for(int m=0;m<nSeason;m++){
                w0+=Math.pow(1-(m+1-0.35)/nSeason,0)*seasonXaccsort[m]/nSeason;
                w1+=Math.pow(1-(m+1-0.35)/nSeason,1)*seasonXaccsort[m]/nSeason;
                w2+=Math.pow(1-(m+1-0.35)/nSeason,2)*seasonXaccsort[m]/nSeason;
            }
            double beta = (2 * w1 - w0)/(6 * w1-w0-6 * w2);
            double g1 = Gamma.gamma(1+1.0/beta);
            double g2 = Gamma.gamma(1-1.0/beta);
            double alpha = (w0-2*w1)*beta/(g1*g2);
            double litteGamma = w0-alpha*g1*g2;
            for (int i=0;i<yearnum;i++){
                double fx = Math.pow(1+Math.pow(alpha/(xacc[i][j]-litteGamma),beta),-1);
                double P = 1-fx;
                double W = 0;
                double c0 = 2.515517;
                double c1 = 0.802853;
                double c2 = 0.010328;
                double d1 = 1.432788;
                double d2 = 0.189269;
                double d3 = 0.001308;
                if(P <= 0.5){
                    W = Math.sqrt(-2.0*Math.log(P));
                    spei[i][j] = W-(c0 + c1 *W+ c2 *W*W)/(1+ d1 *W+ d2 *W*W+ d3 *W*W*W);
                }
                else{
                    W = Math.sqrt(-2.0*Math.log(1-P));
                    spei[i][j] = (c0 + c1 *W+ c2 *W*W)/(1+ d1 *W+ d2 *W*W+ d3 *W*W*W)-W;
                }
            }
        }
        return spei;
    }

    private double[] bubbleSort(double[] arr) {
        for(int i = 0; i < arr.length-1; i++) {
            for(int j = 0; j < arr.length - 1 - i; j++) {
                if(arr[j] > arr[j + 1]){
                    double temp = arr[j];
                    arr[j] = arr[j + 1];
                    arr[j + 1] = temp;
                }
            }
        }
        return arr;
    }

    private double[][] calculateX(){
        double[][] xacc=new double[yearnum][monthnum];
        if (k==1){
            xacc=mothlyD;
        }
        else if(k<=12){
            double[] seriesD = new double[yearnum*monthnum];
            int kk=0;
            for(int i = 0;i<yearnum;i++){
                for(int j=0;j<monthnum;j++){
                    seriesD[kk]=mothlyD[i][j];
                    kk++;
                }
            }
            int seriesLength = kk;
            double[] seriesXacc = new double[yearnum*monthnum];
            for(int i=0;i<seriesLength;i++){
                int backstep = 0;
                while(i-backstep>=0){
                    seriesXacc[i] += seriesD[i-backstep];
                    backstep++;
                    if (backstep==this.k)
                        break;
                }
            }

            for(int i = 0;i<yearnum;i++){
                for(int j=0;j<monthnum;j++){
                    xacc[i][j] = seriesXacc[i*12+j];
                }
            }
        }
        return xacc;
    }
}

class YearMonthlyDI{
    /*
     * 计算某年 每月的D: 降雨-潜在蒸散发
     * 单点角度
     *
     * */
    private final double[][] temp;
    private final double[][] precipitation;
    private final double latitude;
    private final int startyear;
    private final int yearnum;
    private final int monthnum;

    public YearMonthlyDI(double[][] temp,double[][] precipitation,double latitude,int startyear){
        // 时间维度1
        /*
         *  temp 温度数据
         *  precipitation 降水数据
         *  latitude 纬度
         *
         * */
        this.temp = temp;
        this.precipitation = precipitation;
        this.latitude = latitude;
        this.startyear = startyear;
        this.yearnum = temp.length;
        this.monthnum = temp[0].length;
    }

    public double[][] getResult(){
        return calculateD();
    }
    // 计算 D
    private double[][] calculateD(){
        /*************计算PET*************/

        // heat index 为12个月的i的和
        double index_I = 0;  //1
        for(int i=0;i<monthnum;i++){
            double T=0,k=0;
            for(int j=0;j<yearnum;j++)
            {
                T+=temp[j][i];
                k++;
            }
            T=T/k;
            if(T>0)
                index_I += headIndex_i(T);
        }

        // m是关于I的函数
        double m;
        m = 6.75e-7*Math.pow(index_I,3)-7.71e-5*index_I*index_I+0.01792*index_I+0.49239;

        double delta;
        double omegas;
        int NDM;
        double[][] K = new double[yearnum][monthnum];   // 12 correction cofficient for PET

        for(int year=0;year<yearnum;year++){
            for(int month=0;month<monthnum;month++){
                int month2 = month+1;
                switch (month2){
                    case 1:
                    case 3:
                    case 5:
                    case 7:
                    case 8:
                    case 10:
                    case 12:
                    {
                        NDM = 31;
                        double J=0;
                        for(int day=1;day<=31;day++){
                            J+=date2DOY(year+startyear,month2,day,10);
                        }
                        J=1.0*J/31;
                        delta = 0.4093*Math.sin(2*Math.PI*J/365-1.405);
                        omegas=Math.acos(-1.0*Math.tan(latitude)*Math.tan(delta));
                        K[year][month] = (omegas*(24.0/Math.PI)/12)*(1.0*NDM/30);
                        break;
                    }
                    case 4:
                    case 6:
                    case 9:
                    case 11:
                    {
                        NDM = 30;
                        double J=0;
                        for(int day=1;day<=30;day++){
                            J+=date2DOY(year+startyear,month2,day,10);
                        }
                        J=1.0*J/30;
                        delta = 0.4093*Math.sin(2*Math.PI*J/365-1.405);
                        omegas=Math.acos(-1.0*Math.tan(latitude)*Math.tan(delta));
                        K[year][month] = (omegas*(24.0/Math.PI)/12)*(1.0*NDM/30);
                        break;
                    }
                    case 2:
                    {
                        double J=0;

                        if((year % 4 == 0 && year % 100 != 0) || (year%400==0 && year % 3200 != 0) || year % 172800 == 0){
                            //闰年
                            NDM = 29;
                            for(int day=1;day<=29;day++){
                                J+=date2DOY(year+startyear,month2,day,10);
                            }
                            J=1.0*J/29;
                            delta = 0.4093*Math.sin(2*Math.PI*J/365-1.405);
                            omegas=Math.acos(-1.0*Math.tan(latitude)*Math.tan(delta));
                            K[year][month] = (omegas*(24.0/Math.PI)/12)*(1.0*NDM/30);
                        }
                        else{
                            NDM = 29;
                            for(int day=1;day<=28;day++){
                                J+=date2DOY(year+startyear,month2,day,10);
                            }
                            J=1.0*J/28;
                            delta = 0.4093*Math.sin(2*Math.PI*J/365-1.405);
                            omegas=Math.acos(-1.0*Math.tan(latitude)*Math.tan(delta));
                            K[year][month] = (omegas*(24.0/Math.PI)/12)*(1.0*NDM/30);
                        }
                        break;
                    }
                }
            }
        }



        double[][] D = new double[yearnum][monthnum];
//        double[][] PET = new double[yearnum][monthnum];
        for(int year=0;year<yearnum;year++){
            for(int month=0;month<monthnum;month++){

                if (temp[year][month]<=0){
                    D[year][month] = precipitation[year][month];
                    continue;
                }
                D[year][month] = precipitation[year][month]-16*K[year][month]*Math.pow(10*temp[year][month]/index_I,m);
            }
        }

        return D;
        /*************计算PET*************/
    }

    private double headIndex_i(double temp){
        double index_i;
        index_i = Math.pow(temp/5,1.514);
        return index_i;
    }

    private int date2DOY(int year,int month,int day,int hour){
        //年月日转儒略日
        hour = 20;
        int doy = 0;
        double JD0 = (int)(365.25*(year-1))+(int)(30.6001*(1+13))+1+hour/24.0+1720981.5;

        double JD2;
        if (month<=2)
            JD2 = (int)(365.25*(year-1))+(int)(30.6001*(month+13))+day+hour/24.0+1720981.5;
        else
            JD2 =(int)(365.25 * year) + (int)(30.6001 * (month + 1)) + day + hour / 24.0 + 1720981.5;

        doy = (int)(JD2- JD0+1);
        return doy;
    }
}
