import java.io.*;
import java.math.BigDecimal;
import  java.util.*;

public class AntiSensePolyA{

  
    
    public   static void main(String args[])
   {
        String KnownGeneFile=args[0];
        String PlaFile=args[1];
        //String Sample=args[2];
        String OutFile="AntiSense_"+args[1];
        String ExpFile="_Sense-AntiSense_gene.out";

        PrintWriter outputStream1=null;

                       try
                      {
                              outputStream1=new PrintWriter(new FileOutputStream(OutFile,false));

                       }

                        catch(FileNotFoundException e)
                        {
                             System.out.println("Error opening the file stuff.txt.");
                               System.exit(0);
                        }
              PrintWriter outputStream2=null;

                       try
                      {
                              outputStream2=new PrintWriter(new FileOutputStream(ExpFile,false));

                       }

                        catch(FileNotFoundException e)
                        {
                             System.out.println("Error opening the file stuff.txt.");
                               System.exit(0);
                        }

           Map[] chr_fs = new HashMap[300];
           Map[] chr_fe = new HashMap[300];
           Map[] chr_rs = new HashMap[300];
           Map[] chr_re = new HashMap[300];
           ArrayList[] Position_f=new ArrayList[300];
           ArrayList[] Position_r=new ArrayList[300];
           Map Anti_gene=new HashMap();


           for(int i=0;i<300;i++)
           {
               chr_fs[i]=new HashMap();
               chr_fe[i]=new HashMap();
               chr_rs[i]=new HashMap();
               chr_re[i]=new HashMap();
               Position_f[i]=new ArrayList();
               Position_r[i]=new ArrayList();
           }

           Map chr_NO=new HashMap();
           Map NO_chr=new HashMap();
           int chr_no=0;
            try
       {
            BufferedReader inputStream=new BufferedReader(new FileReader(KnownGeneFile));
             String nextline1=new String();
            nextline1=inputStream.readLine();
            nextline1=inputStream.readLine();

            while(nextline1!=null)
            {
                String chr=getTabChar(nextline1,1,2,'	');
                String gene=getTabChar(nextline1,0,1,'	');
                String orientation=getTabChar(nextline1,2,3,'	');
                String txStart=getTabChar(nextline1,3,4,'	');
                String txEnd=getTabChar(nextline1,4,5,'	');
                //String exon_intron=getTabChar(nextline1,8,10,'	');
                String transcript=getTabChar(nextline1,10,11,'	');
                if(!transcript.equals("n/a"))
                {
                     int x;
                    if(chr_NO.containsKey(chr))
                        x=Integer.parseInt(chr_NO.get(chr).toString());
                    else
                      {

                          chr_NO.put(chr, Integer.toString(chr_no));
                          NO_chr.put(Integer.toString(chr_no), chr);
                          x=chr_no;
                          chr_no=chr_no+1;

                       }
                     if(orientation.equals("+"))
                     {
                         chr_fs[x].put(txStart, nextline1);
                         chr_fe[x].put(txEnd, txStart);
                         if(!Position_f[x].contains(txStart)&&!Position_f[x].contains(txEnd))
                         {
                             Position_f[x].add(Integer.parseInt(txStart));
                             Position_f[x].add(Integer.parseInt(txEnd));
                         }
                        
                     }
                     else
                     {
                         chr_rs[x].put(txStart, nextline1);
                         chr_re[x].put(txEnd, txStart);
                         if(!Position_r[x].contains(txStart)&&!Position_r[x].contains(txEnd))
                         {
                             Position_r[x].add(Integer.parseInt(txStart));
                             Position_r[x].add(Integer.parseInt(txEnd));
                         }

                     }
                }

                 nextline1=inputStream.readLine();

            }
            inputStream.close();

        }

              catch(FileNotFoundException e)
         {
             System.out.println("File not found");
             System.out.println("or could not be opened");

         }
       catch (IOException e)
       {
           System.out.println(e);
           System.out.println("Error reading !");
       }
           for(int i=0;i<chr_NO.size();i++)
           {
               Object tmp[]=Position_f[i].toArray();
               Arrays.sort(tmp);
//               for(int x=0;x<10;x++)
//                   System.out.println(tmp[x]);
                for(int x=0;x<tmp.length;x++)
                    Position_f[i].set(x, tmp[x]);
               
               tmp=Position_r[i].toArray();
                Arrays.sort(tmp);
                for(int x=0;x<tmp.length;x++)
                    Position_r[i].set(x, tmp[x]);

                //analysis known sense-antisense
                for(int x=0;x<Position_r[i].size();x++)
                {
                    int search=Integer.parseInt((Position_r[i].get(x).toString()));
                    //System.out.println(search);
                    List ss;
                    String s1=new String();
                    String s2=new String();
                    if(Position_f[i].size()>0)
                    {
                       ss=binary_search(search,Position_f[i]);
                       s1=ss.get(0).toString();
                       s2=ss.get(1).toString();
                    }
                    //s1<=search<=s2
                    //If s1,s2 are the txStart,txEnd of the same gene, then the two genes are sense-antisense genes
                    if(chr_fs[i].containsKey(s1)&&chr_fe[i].containsKey(s2)&&chr_fe[i].get(s2).toString().equals(s1))
                    {
                        String line1=chr_fs[i].get(s1).toString();
                        String line2=new String();
                        if(chr_rs[i].containsKey(Integer.toString(search)))
                            line2=chr_rs[i].get(Integer.toString(search)).toString();
                        else if(chr_re[i].containsKey(Integer.toString(search)))
                            line2=chr_rs[i].get(chr_re[i].get(Integer.toString(search))).toString();
                        String gene1=getTabChar(line1,0,1,'	');
                        String gene2=getTabChar(line2,0,1,'	');
                        if(!Anti_gene.containsKey(gene1))
                            Anti_gene.put(gene2, s1);

                    }
                    //If s1,s2 are the txEnd,txStart of the different genes, then no action
                }



           }

                try
       {
            BufferedReader inputStream=new BufferedReader(new FileReader(PlaFile));
             String nextline1=new String();
            nextline1=inputStream.readLine();
            nextline1=inputStream.readLine();

            while(nextline1!=null)
            {
                String chr=getTabChar(nextline1,1,2,'	');
                String orientation=getTabChar(nextline1,2,3,'	');
                String gene=getTabChar(nextline1,4,5,'	');
                int site=Integer.parseInt(getTabChar(nextline1,3,4,'	'));
                if(orientation.equals("-"))
                {
                    if(chr_NO.containsKey(chr))
                    {
                        int i=Integer.parseInt(chr_NO.get(chr).toString());
                        List ss=binary_search(site,Position_f[i]);
                        String s1=ss.get(0).toString();
                        String s2=ss.get(1).toString();
                        if(chr_fs[i].containsKey(s1)&&chr_fe[i].containsKey(s2)&&chr_fe[i].get(s2).toString().equals(s1))
                    {
                        String line1=chr_fs[i].get(s1).toString();
                        String gene1=getTabChar(line1,0,1,'	');
                        String txStart=getTabChar(line1,3,4,'	');
                        String txEnd=getTabChar(line1,4,5,'	');
                        String cdsStart=getTabChar(line1,5,6,'	');
                        String cdsEnd=getTabChar(line1,6,7,'	');
                        String exon_s=getTabChar(line1,8,9,'	');
                        String exon_e=getTabChar(line1,9,10,'	');
                        String locus_type=new String();
                        Map cds_s=new HashMap();
                        Map cds_e=new HashMap();
                        ArrayList cd=new ArrayList();
                        while(exon_s.indexOf(',')!=-1)
                        {
                            String x1=exon_s.substring(0,exon_s.indexOf(','));
                            exon_s=exon_s.substring(exon_s.indexOf(',')+1);
                            if(!x1.equals(txStart))
                                cds_s.put(x1, "0");
                            cd.add(Integer.parseInt(x1));
                        }
                        while(exon_e.indexOf(',')!=-1)
                        {
                            String x1=exon_e.substring(0,exon_e.indexOf(','));
                            exon_e=exon_e.substring(exon_e.indexOf(',')+1);
                            if(!x1.equals(txEnd))
                                cds_e.put(x1, "0");
                            cd.add(Integer.parseInt(x1));
                        }
                        cd.add(Integer.parseInt(cdsStart));
                        cd.add(Integer.parseInt(cdsEnd));

                         Object tmp[]=cd.toArray();
                        Arrays.sort(tmp);
                         for(int x=0;x<tmp.length;x++)
                             cd.set(x, tmp[x]);
                        ss=binary_search(site,cd);
                        s1=ss.get(0).toString();
                        s2=ss.get(1).toString();

                        int number=0;
                        for(int x=0;x<tmp.length;x++)
                        {
                            if(cd.get(x).equals(s2))
                            {
                                number=x;
                                break;

                            }
                        }
                        if(number==1)
                            locus_type="5UTR";
                        else if(number==cd.size()-1)
                            locus_type="3UTR";
                        else if(number%2==0)
                            locus_type="exon";
                        else
                            locus_type="intron";

                        outputStream1.println(nextline1+'	'+gene1+'	'+locus_type);      

                    }
                        else if(chr_fe[i].containsKey(s1)&&chr_fs[i].containsKey(s2))
                        {
                            if(site-Integer.parseInt(s1)>Integer.parseInt(s2)-site)
                            {
                                if(Integer.parseInt(s2)-site<1000)
                                {
                                    String locus_type="upstream1k";
                                    String gene1=getTabChar(chr_fs[i].get(s2).toString(),0,1,'	');
                                     outputStream1.println(nextline1+'	'+gene1+'	'+locus_type);

                                }

                            }
                            else
                            {
                                if(site-Integer.parseInt(s1)<1000)
                                {
                                    String locus_type="downstream1k";
                                    String gene1=getTabChar(chr_fs[i].get(chr_fe[i].get(s1)).toString(),0,1,'	');
                                     outputStream1.println(nextline1+'	'+gene1+'	'+locus_type); 
                                }

                            }
                        }


                    }
                }
                nextline1=inputStream.readLine();

            }
            inputStream.close();

        }

              catch(FileNotFoundException e)
         {
             System.out.println("File not found");
             System.out.println("or could not be opened");

         }
       catch (IOException e)
       {
           System.out.println(e);
           System.out.println("Error reading !");
       }
           
           
           
           outputStream1.close();
           outputStream2.close();






     
    }

    public static String getTabChar(String s,int m,int n,char label)
    {

        int block=0;
        int locus1=-1,locus2=s.length();
        for(int x=0;x<s.length();x++)
        {
            if(s.charAt(x)==label)
            {
                  block=block+1;
                 if(block==m) locus1=x;
                 else if(block==n)  {locus2=x;break;}
            }
         }
        return s.substring(locus1+1, locus2);
    }

    public static List binary_search(int x, ArrayList y)
    {
        List t=y;
        while(t.size()>2)
        {
            if(Integer.parseInt(t.get(t.size()/2).toString())>x)
                t=t.subList(0,t.size()/2+1);
            else
                t=t.subList(t.size()/2, t.size());
        }
          return t;
    }
    
}

