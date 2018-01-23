/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package activecontourtracker;

import activecontourtracker.ContourTracker.Contour;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonSyntaxException;
import java.awt.AWTException;
import java.awt.AlphaComposite;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.awt.geom.GeneralPath;
import java.awt.image.AffineTransformOp;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferByte;
import java.awt.image.DataBufferUShort;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.DebugTools;
import loci.formats.FormatException;
import loci.formats.ImageReader;
import loci.formats.Memoizer;
import loci.formats.meta.MetadataStore;
import loci.formats.ome.OMEXMLMetadataImpl;

/**
 *
 * @author edwardharry
 */
class Dims
{

    static
            Dims Zero()
    {
        return new Dims(0, 0);
    }

    int width, height;

    Dims(int w, int h)
    {
        width = w;
        height = h;
    }

}

public
        class ActiveContourTracker
{

    ///////////////////////  PARAMS START  ///////////////////////
    private final
            double rollingBallBgRadius;
    private final
            double medianFilterKernelRadius;
    private final
            double erodeFilterKernelRadius;
    private final
            int nErodeFilterPasses;

    private final
            double contourSampling;
    private final
            double contourDivisionSensitivity;
    private final
            double contourConvCritSq;
    private final
            int contourSlidingWindowSize;
    private final
            double contourTimeStep;
    private final
            int contourConvergence_nIter;

    private final
            double internalWeight;
    private final
            double edgeWeight;
    private final
            double regionWeight;
    private final
            double regionSensitivity;
    private final
            double axisWeight;
    private final
            double balloonWeight;
    private final
            double volumeConstraint;
    private final
            boolean couplingFlag;

    private final
            double maskMovieImageScale;
    private final
            double summaryMovieImageScale;
    private final
            int movieFrameRate;
    private final
            int movieSimpleFrameEffectiveFrameRate;
    private final
            int movieDetailedFrameInterval;
    ///////////////////////  PARAMS  END   ///////////////////////

    private final
            ChannelLookup channelLookup;
    private final
            Processor multithreader;
    private final
            ImageReaderPool imageReaderPool;
    private final
            Dims dims;
    private final
            RollingBallBackgroundSubtractor.RollingBall ball;
    private final
            BinaryMask.FilterKernel medianFilterKernel;
    private final
            BinaryMask.FilterKernel erodeFilterKernel;
    private final
            BinaryMask.WatershedParams watershedParams;
    private final
            String nucleiMovieFile;
    private final
            String wholeCellMovieFile;
    private final
            String summaryMovieFile;
    private final
            int seriesIndex;
    private final
            int nTimePoints;
    private final
            HDF5Writer h5Writer;
    private final
            ArrayList<ROIOutput> roiOutputs;
    private final
            ArrayList<ChOutput> chOutputs;

    private final static
            String[] CHANNEL_OUTPUT_ORDER =
            {
                "mCherry-FUCCI", "CFP 425-FUCCI",
                "Hyper 425-FUCCI", "Hyper 500-fucci"
            };

    private final static
            String SETTINGS_FILE_NAME = "ActiveContourSettings.json";

    private static
            class SettingsLoader
    {

        static
                class Settings
        {

            public
                    double rollingBallBgRadius;
            public
                    double medianFilterKernelRadius;
            public
                    double erodeFilterKernelRadius;
            public
                    int nErodeFilterPasses;

            public
                    double contourSampling;
            public
                    double contourDivisionSensitivity;
            public
                    double contourConvCritSq;
            public
                    int contourSlidingWindowSize;
            public
                    double contourTimeStep;
            public
                    int contourConvergence_nIter;

            public
                    double internalWeight;
            public
                    double edgeWeight;
            public
                    double regionWeight;
            public
                    double regionSensitivity;
            public
                    double axisWeight;
            public
                    double balloonWeight;
            public
                    double volumeConstraint;
            public
                    boolean couplingFlag;

            public
                    boolean recordMaskMovies;
            public
                    boolean recordSummaryMovie;
            public
                    double maskMovieImageScale;
            public
                    double summaryMovieImageScale;
            public
                    int movieFrameRate;
            public
                    int movieSimpleFrameEffectiveFrameRate;
            public
                    int movieDetailedFrameInterval;
        }

        static
                class DefaultSettings extends Settings
        {

            DefaultSettings()
            {
                movieDetailedFrameInterval = 20;
                movieSimpleFrameEffectiveFrameRate = 10;
                movieFrameRate = 30;
                summaryMovieImageScale = 0.35;
                maskMovieImageScale = 0.7;
                couplingFlag = true;
                volumeConstraint = 0.01;
                balloonWeight = 0.01;
                axisWeight = 0;
                regionSensitivity = 0.1;
                regionWeight = 0.001;
                edgeWeight = 0.1;
                internalWeight = 0.01;
                contourConvergence_nIter = 100;
                contourTimeStep = 0.5;
                contourSlidingWindowSize = 20;
                contourConvCritSq = 0.01;
                contourDivisionSensitivity = 0.1;
                contourSampling = 4.0;
                nErodeFilterPasses = 3;
                erodeFilterKernelRadius = 1.0;
                medianFilterKernelRadius = 5.0;
                rollingBallBgRadius = 150.0;
                recordMaskMovies = false;
                recordSummaryMovie = false;
            }
        }

        static
                Settings Load(File file)
        {
            Settings settings = null;

            if (file.canRead())
            {
                try
                {
                    settings = new GsonBuilder().create().fromJson(new FileReader(file), Settings.class);
                    System.out.println("Read settings from file " + file.getPath());
                }
                catch (FileNotFoundException ex)
                {
                    Logger.getLogger(ActiveContourTracker.class.getName()).
                            log(Level.SEVERE, null, ex);
                }
                catch (JsonSyntaxException ex)
                {
                    System.err.println(ex.getMessage());
                    System.err.println("Error with settings file " + file.getPath()
                                       + ". Using default parameters");

                    settings = new DefaultSettings();
                }
            }
            else
            {
                System.out.println("Using default parameters");
                settings = new DefaultSettings();
            }

            return settings;
        }

    }

    /**
     * @param args the command line arguments
     */
    public static
            void main(String[] args)
    {
        try
        {
            ActiveContourTracker tracker = new ActiveContourTracker(args[0],
                                                                    Integer.parseUnsignedInt(args[1]));
            tracker.ReadSeries();
        }
        catch (FormatException | IOException | InterruptedException | ExecutionException ex)
        {
            Logger.getLogger(ActiveContourTracker.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public
            ActiveContourTracker(String fileString, int seriesIndex)
            throws IOException
    {
        DebugTools.enableLogging("INFO");

        File file = new File(fileString);

        SettingsLoader.Settings settings = SettingsLoader.Load(new File(file.getParent(),
                                                                        SETTINGS_FILE_NAME));

        rollingBallBgRadius = settings.rollingBallBgRadius;
        medianFilterKernelRadius = settings.medianFilterKernelRadius;
        erodeFilterKernelRadius = settings.erodeFilterKernelRadius;
        nErodeFilterPasses = settings.nErodeFilterPasses;
        contourSampling = settings.contourSampling;
        contourDivisionSensitivity = settings.contourDivisionSensitivity;
        contourConvCritSq = settings.contourConvCritSq;
        contourSlidingWindowSize = settings.contourSlidingWindowSize;
        contourTimeStep = settings.contourTimeStep;
        contourConvergence_nIter = settings.contourConvergence_nIter;
        internalWeight = settings.internalWeight;
        edgeWeight = settings.edgeWeight;
        regionWeight = settings.regionWeight;
        regionSensitivity = settings.regionSensitivity;
        axisWeight = settings.axisWeight;
        balloonWeight = settings.balloonWeight;
        volumeConstraint = settings.volumeConstraint;
        couplingFlag = settings.couplingFlag;
        maskMovieImageScale = settings.maskMovieImageScale;
        summaryMovieImageScale = settings.summaryMovieImageScale;
        movieFrameRate = settings.movieFrameRate;
        movieSimpleFrameEffectiveFrameRate = settings.movieSimpleFrameEffectiveFrameRate;
        movieDetailedFrameInterval = settings.movieDetailedFrameInterval;

        ball = new RollingBallBackgroundSubtractor.RollingBall(rollingBallBgRadius);
        medianFilterKernel = new BinaryMask.FilterKernel(medianFilterKernelRadius);
        erodeFilterKernel = new BinaryMask.FilterKernel(erodeFilterKernelRadius);

        multithreader = new Processor();

        System.out.println("Started processor with "
                           + multithreader.getCorePoolSize() + " threads");

        imageReaderPool = new ImageReaderPool(fileString, seriesIndex);

        Memoizer reader = imageReaderPool.GetImageReader();
        channelLookup = new ChannelLookup(reader.getMetadataStore());
        dims = new Dims(reader.getSizeX(), reader.getSizeY());
        nTimePoints = reader.getSizeT();
        imageReaderPool.AddImageReader(reader);

        System.out.println("" + nTimePoints + " timepoints to "
                           + "process for series " + seriesIndex);

        watershedParams = new BinaryMask.WatershedParams(dims.width);

        String name = file.getName();
        if (name != null && name.contains("."))
        {
            name = name.substring(0, name.lastIndexOf('.'));
        }

        h5Writer = new HDF5Writer(new File(file.getParent(), name + "_XY"
                                                             + seriesIndex + ".h5").getPath());

        if (settings.recordMaskMovies)
        {
            nucleiMovieFile = new File(file.getParent(), name + "_XY"
                                                         + seriesIndex
                                                         + "_Nuclei.webm").getPath();
            wholeCellMovieFile = new File(file.getParent(), name + "_XY"
                                                            + seriesIndex
                                                            + "_WholeCell.webm").getPath();
        }
        else
        {
            nucleiMovieFile = wholeCellMovieFile = null;
        }

        if (settings.recordSummaryMovie)
        {
            summaryMovieFile = new File(file.getParent(), name + "_XY"
                                                          + seriesIndex
                                                          + "_Summary.webm").getPath();
        }
        else
        {
            summaryMovieFile = null;
        }

        this.seriesIndex = seriesIndex;

        roiOutputs = new ArrayList<>();
        chOutputs = new ArrayList<>();
        CreateOutputs();
        PrepareOutputDataFunctions();
    }

    private
            void CreateOutputs()
    {
        AddROIOutput("Centre_X", ((roi) ->
        {
            return roi.centre.x;
        }));

        AddROIOutput("Centre_Y", ((roi) ->
        {
            return roi.centre.y;
        }));

        AddROIOutput("Perimeter", ((roi) ->
        {
            return roi.perimeter;
        }));

        AddROIOutput("Area", ((roi) ->
        {
            return roi.area;
        }));

        AddROIOutput("Track_ID", ((roi) ->
        {
            return roi.trackID;
        }));

        AddROIOutput("Previous_Track_ID", ((roi) ->
        {
            return roi.previousTrackID;
        }));

        AddChOutput("Total_Intensity", ((stat) ->
        {
            return stat.total;
        }));

        AddChOutput("Mean_Intensity", ((stat) ->
        {
            return stat.mean;
        }));

        AddChOutput("SD_Intensity", ((stat) ->
        {
            return stat.sd;
        }));
    }

    @FunctionalInterface
    interface ROIfunction<R>
    {

        R apply(ROIExtractor.ROI roi);
    }

    @FunctionalInterface
    interface Statfunction<R>
    {

        R apply(ROIStatistics.Statistics stat);
    }

    private
            class ROIOutput
    {

        final
                String name;
        final
                ROIfunction func;

        ROIOutput(String name, ROIfunction func)
        {
            this.name = name;
            this.func = func;
        }
    }

    private
            class ChOutput
    {

        final
                String name;
        final
                Statfunction func;

        ChOutput(String name, Statfunction func)
        {
            this.name = name;
            this.func = func;
        }
    }

    private
            void AddROIOutput(String name, ROIfunction func)
    {
        roiOutputs.add(new ROIOutput(name, func));
    }

    private
            void AddChOutput(String name, Statfunction func)
    {
        chOutputs.add(new ChOutput(name, func));
    }

    private
            enum ROIType
    {
        Nucleus("Nucleus"),
        Cell("WholeCell");

        private final
                String str;

        ROIType(String str)
        {
            this.str = str;
        }

        @Override
        public
                String toString()
        {
            return str;
        }

        static
                String GetLongestName()
        {
            String out = "";
            for (ROIType typ : ROIType.values())
            {
                String typS = typ.toString();

                if (typS.length() > out.length())
                {
                    out = typS;
                }
            }

            return out;
        }
    }

    private
            String[] outputNames;
    private
            Object[] outputTypes;
    private
            boolean outputTypesSet;
    private
            ROIfunction[] roiFuncs;
    private
            Statfunction[] statFuncs;

    private
            Object GetTypeRep(Object in)
    {
        Object out = null;
        if (in instanceof Number)
        {
            if (in instanceof Double)
            {
                out = Double.NaN;
            }
            else if (in instanceof Float)
            {
                out = Float.NaN;
            }
            else if (in instanceof Long)
            {
                out = Long.MAX_VALUE;
            }
            else if (in instanceof Short)
            {
                out = Short.MAX_VALUE;
            }
            else if (in instanceof Integer)
            {
                out = Integer.MAX_VALUE;
            }
        }
        else if (in instanceof String)
        {
            out = in;
        }

        if (out == null)
        {
            throw new RuntimeException("Type " + in + " not supported for data output");
        }

        return out;
    }

    private
            int PrepareOutputHeader(boolean fillArrays)
    {
        int index = 0;

        if (fillArrays)
        {
            outputNames[index] = "SeriesIndex";
            outputTypes[index] = GetTypeRep(0);
        }
        ++index;

        if (fillArrays)
        {
            outputNames[index] = "TimePoint";
            outputTypes[index] = GetTypeRep(0);
        }
        ++index;

        if (fillArrays)
        {
            outputNames[index] = "ROI_Type";
            outputTypes[index] = GetTypeRep(ROIType.GetLongestName());
        }
        ++index;

        return index;
    }

    private
            int CreateOutputHeader(Object[] col, final
                                   int series,
                                   final
                                   int time, final
                                   String name)
    {
        int index = 0;

        col[index++] = series;
        col[index++] = time;
        col[index++] = name;

        return index;
    }

    private
            void PrepareOutputDataFunctions()
    {
        outputTypesSet = false;

        final
                int nBase = PrepareOutputHeader(false);
        final
                int nCh = CHANNEL_OUTPUT_ORDER.length;
        final
                int nROI = roiOutputs.size();
        final
                int nPerCh = chOutputs.size();

        final
                int nTotal = nBase + nROI + (nPerCh * nCh);

        outputNames = new String[nTotal];
        outputTypes = new Object[nTotal];
        roiFuncs = new ROIfunction[nROI];
        statFuncs = new Statfunction[nPerCh];

        PrepareOutputHeader(true);

        int index = nBase;

        for (int roiIndex = 0;
             roiIndex < nROI;
             ++roiIndex, ++index)
        {
            outputNames[index] = roiOutputs.get(roiIndex).name;
            roiFuncs[roiIndex] = roiOutputs.get(roiIndex).func;
        }

        for (int chIndex = 0;
             chIndex < nCh;
             ++chIndex)
        {
            for (int statIndex = 0;
                 statIndex < nPerCh;
                 ++statIndex, ++index)
            {
                outputNames[index] = chOutputs.get(statIndex).name + "_"
                                     + CHANNEL_OUTPUT_ORDER[chIndex];

                if (chIndex == 0)
                {
                    statFuncs[statIndex] = chOutputs.get(statIndex).func;
                }
            }
        }

        roiOutputs.clear();
        chOutputs.clear();
    }

    private
            Object[][] CreateOutputData(ROIExtractor.ROI[] rois, ROIType type,
                                        ROIStatistics stats)
    {
        final
                int nCol = outputNames.length;
        final
                int nRow = rois.length;
        final
                int nCh = CHANNEL_OUTPUT_ORDER.length;

        ArrayList<Object[]> data = new ArrayList<>(nRow);

        for (int iRow = 0;
             iRow < nRow;
             ++iRow)
        {
            ROIExtractor.ROI roi = rois[iRow];

            if (roi == null || stats.stats[0][iRow] == null)
            {
                continue;
            }

            Object[] col = new Object[nCol];

            int index = CreateOutputHeader(col, seriesIndex,
                                           roi.t, type.toString());
            final
                    int nBase = index;

            for (int roiIndex = 0;
                 roiIndex < roiFuncs.length;
                 ++roiIndex)
            {
                col[index++] = roiFuncs[roiIndex].apply(roi);
            }

            for (int chIndex = 0;
                 chIndex < nCh;
                 ++chIndex)
            {
                for (int statIndex = 0;
                     statIndex < statFuncs.length;
                     ++statIndex)
                {
                    col[index++] = statFuncs[statIndex].apply(stats.stats[chIndex][iRow]);
                }
            }

            if (outputTypesSet == false)
            {
                for (int colIndex = nBase;
                     colIndex < col.length;
                     ++colIndex)
                {
                    outputTypes[colIndex] = GetTypeRep(col[colIndex]);
                }

                outputTypesSet = true;
            }

            data.add(col);
        }

        return data.toArray(new Object[0][0]);
    }

    private
            class ImageReaderPool
    {

        private final
                ArrayDeque<Memoizer> imageReaders;

        public
                ImageReaderPool(final
                        String file, final
                                int seriesIndex)
        {
            int n = multithreader.getCorePoolSize();
            imageReaders = new ArrayDeque<>(n);
            Future<?>[] tasks = new Future<?>[n];

            for (int index = 0;
                 index < n;
                 ++index)
            {
                Runnable task = ()
                        ->
                {
                    try
                    {
                        final
                                Memoizer imageReader = new Memoizer(new ImageReader());
                        imageReader.setId(file);
                        imageReader.setSeries(seriesIndex);
                        synchronized (imageReaders)
                        {
                            imageReaders.add(imageReader);
                        }
                    }
                    catch (IOException | FormatException ex)
                    {
                        System.err.println(ex.getMessage());
                    }
                };
                tasks[index] = multithreader.submit(task);
            }

            for (Future<?> future : tasks)
            {
                try
                {
                    future.get();
                }
                catch (InterruptedException | ExecutionException ex)
                {
                    System.err.println(ex.getMessage());
                }
            }
        }

        private
                Memoizer GetImageReader()
        {
            final
                    Memoizer reader;
            synchronized (imageReaders)
            {
                reader = imageReaders.poll();
            }
            return reader;
        }

        private
                void AddImageReader(final
                        Memoizer reader)
        {
            synchronized (imageReaders)
            {
                imageReaders.add(reader);
            }
        }

        private
                void Shutdown()
        {
            imageReaders.stream().forEach((reader)
                    ->
            {
                try
                {
                    reader.close();
                }
                catch (IOException ex)
                {
                    System.err.println(ex.getMessage());
                }
            });
        }

    }

    short[] ReadImage(Memoizer reader, String chName, int timePoint)
            throws FormatException, IOException
    {
        return ReadImage(reader, reader.getIndex(0,
                                                 channelLookup.GetChannelIndex(chName), timePoint));
    }

    short[] ReadImage(Memoizer reader, int index) throws FormatException, IOException
    {
        short[] image = new short[dims.width * dims.height];
        
        switch (reader.getBitsPerPixel())
        {
            case 16:
                ByteBuffer.wrap(reader.openBytes(index)).order(reader.isLittleEndian()
                                                       ? ByteOrder.LITTLE_ENDIAN : ByteOrder.BIG_ENDIAN).
                    asShortBuffer().get(image);
                break;
                
            case 8:
                byte[] im_08 = reader.openBytes(index);
                for (int pixelIndex = 0;
                     pixelIndex < im_08.length;
                     ++pixelIndex)
                {
                    image[pixelIndex] = (short)(((im_08[pixelIndex] & 0xff) / 255.0 * 65535.0) + 0.5);
                }
                break;
        }
        
        return image;
    }

    private
            class Frame
    {

        final
                short[] fucciG1, fucciG2, redox425, redox500;
        final
                byte[] redoxMask, fucciMask, redoxEdgeMask, fucciEdgeMask;

        final
                HashMap<String, short[]> channelSelector;

        public
                Frame(final
                        int t) throws FormatException, IOException,
                                      InterruptedException, ExecutionException
        {
            int nCh = CHANNEL_OUTPUT_ORDER.length;

            ArrayList<Future<short[]>> imageReadTasks = new ArrayList<>(nCh);

            final
                    int readTimeoutMilS = 10000;

            for (int chIndex = 0;
                 chIndex < nCh;
                 ++chIndex)
            {
                final
                        int myId = chIndex;

                imageReadTasks.add(chIndex, multithreader.submit(()
                                   ->
                           {
                               int timer = 0;
                               Memoizer reader = null;
                               RuntimeException exception = null;
                               short[] output = null;

                               try
                               {
                                   while (reader == null && timer < readTimeoutMilS)
                                   {
                                       reader = imageReaderPool.GetImageReader();

                                       if (reader != null)
                                       {
                                           output = ReadImage(reader, CHANNEL_OUTPUT_ORDER[myId], t);
                                           RollingBallBackgroundSubtractor.
                                                   SubtractBackground(output, dims, ball);
                                       }
                                       else
                                       {
                                           Thread.sleep(1);
                                           ++timer;
                                       }
                                   }

                                   if (reader == null)
                                   {
                                       throw new RuntimeException("Timed out waiting"
                                                                  + "for an image reader to become"
                                                                  + " available");
                                   }
                               }
                               catch (FormatException | IOException ex)
                               {
                                   exception = new RuntimeException(ex.getMessage(),
                                                                    ex.getCause());
                               }
                               finally
                               {
                                   imageReaderPool.AddImageReader(reader);
                               }

                               if (exception != null)
                               {
                                   throw exception;
                               }

                               return output;
                           }
                           ));
            }

            List<String> chNames = Arrays.asList(CHANNEL_OUTPUT_ORDER);
            fucciG1 = imageReadTasks.get(chNames.indexOf("mCherry-FUCCI")).get();
            fucciG2 = imageReadTasks.get(chNames.indexOf("CFP 425-FUCCI")).get();
            redox425 = imageReadTasks.get(chNames.indexOf("Hyper 425-FUCCI")).get();
            redox500 = imageReadTasks.get(chNames.indexOf("Hyper 500-fucci")).get();
            imageReadTasks.clear();

            channelSelector = new HashMap<>(nCh);
            channelSelector.put("mCherry-FUCCI", fucciG1);
            channelSelector.put("CFP 425-FUCCI", fucciG2);
            channelSelector.put("Hyper 425-FUCCI", redox425);
            channelSelector.put("Hyper 500-fucci", redox500);

            Future<byte[]> redoxMaskJob = multithreader.submit(()
                    ->
            {
                return BinaryMask.MakeMask(redox500);
            });
            Future<byte[]> fucciG1MaskJob = multithreader.submit(()
                    ->
            {
                return BinaryMask.MakeMask(fucciG1);
            });
            Future<byte[]> fucciG2MaskJob = multithreader.submit(()
                    ->
            {
                return BinaryMask.MakeMask(fucciG2);
            });

            redoxMask = redoxMaskJob.get();
            byte[] fucciG1Mask = fucciG1MaskJob.get();
            byte[] fucciG2Mask = fucciG2MaskJob.get();

            byte[] fucciNuclei = BinaryMask.AddMasks(fucciG1Mask, fucciG2Mask);
            byte[] fucciNoRedox = BinaryMask.SubtractMasks(fucciNuclei, redoxMask);
            fucciMask = BinaryMask.SubtractMasks(fucciNuclei, fucciNoRedox);

            Future<byte[]> fucciMaskFilterJob = multithreader.submit(()
                    ->
            {
                BinaryMask.MedianFilter(fucciMask, dims, medianFilterKernel);
                BinaryMask.Watershed(fucciMask, dims, watershedParams);

                byte[] output = BinaryMask.CopyMask(fucciMask);
                BinaryMask.SobelFilter(output, dims);

                return output;
            });

            Future<byte[]> redoxMaskFilterJob = multithreader.submit(()
                    ->
            {
                BinaryMask.MedianFilter(redoxMask, dims, medianFilterKernel);

                byte[] output = BinaryMask.CopyMask(redoxMask);
                BinaryMask.SobelFilter(output, dims);

                return output;
            });

            fucciEdgeMask = fucciMaskFilterJob.get();
            redoxEdgeMask = redoxMaskFilterJob.get();
        }

        ROIExtractor.ROI[] MakeROIs()
        {
            byte[] fucciROI = BinaryMask.CopyMask(fucciMask);
            for (int index = 0;
                 index < nErodeFilterPasses;
                 ++index)
            {
                BinaryMask.ErodeFilter(fucciROI, dims, erodeFilterKernel);
            }

            return ROIExtractor.SingleImageROIs(fucciROI, dims);
        }
    }

    void CreateTrackers(Frame frame, ContourTracker[] trackers,
                        ContourTracker.ContourUpdate[] updaters)
    {
        ROIExtractor.ROI[] rois = frame.MakeROIs();

        ArrayList<Contour> contour_list = new ArrayList<>(rois.length);
        Future<?>[] tasks = new Future<?>[rois.length];
        int index = -1;
        for (ROIExtractor.ROI roi : rois)
        {
            Runnable task = ()
                    ->
            {
                Contour contour = new Contour(contourSampling,
                                              contourDivisionSensitivity,
                                              contourSlidingWindowSize,
                                              contourConvCritSq, roi);
                synchronized (contour_list)
                {
                    contour_list.add(contour);
                }
            };

            tasks[++index] = multithreader.submit(task);
        }

        for (Future<?> future : tasks)
        {
            try
            {
                future.get();
            }
            catch (InterruptedException | ExecutionException ex)
            {
                System.err.println(ex.getMessage());
            }
        }

        Contour[] contours = contour_list.toArray(new Contour[0]);

        trackers[0] = new ContourTracker(multithreader, contours, dims,
                                         updaters[0], contourSampling, contourDivisionSensitivity,
                                         contourConvCritSq, contourSlidingWindowSize, contourTimeStep,
                                         contourConvergence_nIter, internalWeight, edgeWeight,
                                         regionWeight, regionSensitivity, axisWeight, balloonWeight,
                                         volumeConstraint, couplingFlag, erodeFilterKernel, nErodeFilterPasses);

        trackers[1] = new ContourTracker(multithreader, contours, dims,
                                         updaters[1], contourSampling, contourDivisionSensitivity,
                                         contourConvCritSq, contourSlidingWindowSize, contourTimeStep,
                                         contourConvergence_nIter, internalWeight, edgeWeight,
                                         regionWeight, regionSensitivity, axisWeight, balloonWeight,
                                         volumeConstraint, couplingFlag, erodeFilterKernel, nErodeFilterPasses);
    }

    private
            Future<Frame> FetchFrame(final
                    int t)
    {
        return multithreader.submit(()
                ->
        {
            return new Frame(t);
        });
    }

    void ReadSeries() throws FormatException, IOException,
                             InterruptedException, ExecutionException
    {
        final
                String dataName = "data";

        Future<Frame> nextFrameJob = null;

        Frame frame = new Frame(0);

        if (nTimePoints > 1)
        {
            nextFrameJob = FetchFrame(1);
        }

        MovieFrameGenerator nucleiMovieGen = null;
        MovieFrameGenerator wholeCellMovieGen = null;
        int simpleFrameMultiplier = 0;
        SimpleMovieFrameGenerator movieGen = null;

        if (nucleiMovieFile != null)
        {
            nucleiMovieGen = new MovieFrameGenerator(frame.fucciMask, dims,
                                                     maskMovieImageScale, "Nuclei",
                                                     movieDetailedFrameInterval);

            wholeCellMovieGen = new MovieFrameGenerator(frame.redoxMask, dims,
                                                        maskMovieImageScale, "Whole Cell",
                                                        movieDetailedFrameInterval);
            try
            {
                MovieCreator.Movie nucleiMovie = MovieCreator.VP9Movie(nucleiMovieFile, movieFrameRate,
                                                                       nucleiMovieGen.GetOutputDims());
                nucleiMovieGen.SetUpdater(nucleiMovie);

                MovieCreator.Movie wholeCellMovie = MovieCreator.VP9Movie(wholeCellMovieFile, movieFrameRate,
                                                                          wholeCellMovieGen.GetOutputDims());
                wholeCellMovieGen.SetUpdater(wholeCellMovie);

                simpleFrameMultiplier = Math.max(movieFrameRate
                                                 / movieSimpleFrameEffectiveFrameRate, 1);
            }
            catch (AWTException ex)
            {
                Logger.getLogger(ActiveContourTracker.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        if (summaryMovieFile != null)
        {
            movieGen = new SimpleMovieFrameGenerator(dims, summaryMovieImageScale);
            MovieCreator.Movie movie;
            try
            {
                movie = MovieCreator.VP9Movie(summaryMovieFile, movieSimpleFrameEffectiveFrameRate,
                                              movieGen.GetMovieFrameDims());
                movieGen.SetUpdater(movie);
            }
            catch (AWTException ex)
            {
                Logger.getLogger(ActiveContourTracker.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        ContourTracker.ContourMultiUpdate nucleiUpdater = null;
        ContourTracker.ContourMultiUpdate wholeCellUpdater = null;

        if (nucleiMovieGen != null || movieGen != null)
        {
            ArrayList<ContourTracker.ContourUpdate> nucleiUpdaters
                                                    = new ArrayList<>(2);
            ArrayList<ContourTracker.ContourUpdate> wholeCellUpdaters
                                                    = new ArrayList<>(2);

            if (nucleiMovieGen != null)
            {
                nucleiUpdaters.add(nucleiMovieGen);
            }
            if (movieGen != null)
            {
                nucleiUpdaters.add(movieGen.GetNucleiContours());
            }

            if (wholeCellMovieGen != null)
            {
                wholeCellUpdaters.add(wholeCellMovieGen);
            }
            if (movieGen != null)
            {
                wholeCellUpdaters.add(movieGen.GetWholeCellContours());
            }

            nucleiUpdater = new ContourTracker.ContourMultiUpdate(nucleiUpdaters.toArray(new ContourTracker.ContourUpdate[0]));

            wholeCellUpdater = new ContourTracker.ContourMultiUpdate(wholeCellUpdaters.toArray(new ContourTracker.ContourUpdate[0]));
        }

        ContourTracker fucciTracker;
        ContourTracker redoxTracker;
        {
            ContourTracker[] trackers = new ContourTracker[2];
            CreateTrackers(frame, trackers, new ContourTracker.ContourUpdate[]
                   {
                       nucleiUpdater,
                       wholeCellUpdater
            });
            fucciTracker = trackers[0];
            redoxTracker = trackers[1];
        }

        for (int t = 0;
             t < nTimePoints;
             ++t)

        {
            fucciTracker.EvolveContoursToNextFrame(frame.fucciMask,
                                                   frame.fucciEdgeMask, "Nuclei tracker: ");

            redoxTracker.EvolveContoursToNextFrame(frame.redoxMask,
                                                   frame.redoxEdgeMask, "Whole cell tracker: ");

            {
                ROIExtractor.ROI[] rois = fucciTracker.GetROIs();

                ROIStatistics stats = new ROIStatistics(rois, frame.fucciMask,
                                                        frame);

                h5Writer.WriteCompoundData(dataName, outputNames, outputTypes,
                                           CreateOutputData(rois, ROIType.Nucleus, stats), t == 0);

                rois = redoxTracker.GetROIs();

                stats = new ROIStatistics(rois, frame.redoxMask,
                                          frame);

                h5Writer.WriteCompoundData(dataName, outputNames, outputTypes,
                                           CreateOutputData(rois, ROIType.Cell, stats));
            }

            if (nucleiMovieGen != null)
            {
                for (int index = 0;
                     index < simpleFrameMultiplier;
                     ++index)
                {
                    nucleiMovieGen.PushUpdate();
                }
            }

            if (wholeCellMovieGen != null)
            {
                for (int index = 0;
                     index < simpleFrameMultiplier;
                     ++index)
                {
                    wholeCellMovieGen.PushUpdate();
                }
            }

            if (movieGen != null)
            {
                movieGen.SetFrame(frame);
            }

            if (nextFrameJob != null)
            {
                frame = nextFrameJob.get();

                if (nucleiMovieGen != null)
                {
                    nucleiMovieGen.ChangeImage(frame.fucciMask);
                }

                if (wholeCellMovieGen != null)
                {
                    wholeCellMovieGen.ChangeImage(frame.redoxMask);
                }

                if (t < (nTimePoints - 2))
                {
                    nextFrameJob = FetchFrame(t + 2);
                }
                else
                {
                    nextFrameJob = null;
                }
            }
        }

        System.out.println("Finished, shutting down...");

        if (nucleiMovieGen != null)
        {
            nucleiMovieGen.Close();
        }

        if (wholeCellMovieGen != null)
        {
            wholeCellMovieGen.Close();
        }

        if (movieGen != null)
        {
            movieGen.Close();
        }

        h5Writer.Close();

        imageReaderPool.Shutdown();
        multithreader.ShutdownAndWait();

        System.out.println("Job completed");
    }

    private
            class ROIStatistics
    {

        class Statistics
        {

            String channelName;
            double total;
            double mean;
            double sd;
        }

        final
                Statistics[][] stats;

        public
                ROIStatistics(final
                        ROIExtractor.ROI[] rois, final
                              byte[] mask,
                              final
                              Frame frame)
        {
            final
                    int nChannels = CHANNEL_OUTPUT_ORDER.length;

            stats = new Statistics[nChannels][rois.length];

            Future<?>[] tasks = new Future<?>[rois.length];
            int taskIndex = 0;
            for (ROIExtractor.ROI roi : rois)
            {
                final
                        int myTaskIndex = taskIndex;

                Runnable task = ()
                        ->
                {
                    int nPix = 0;
                    double[] total = new double[nChannels];
                    Arrays.fill(total, 0);
                    double[] total2 = new double[nChannels];
                    Arrays.fill(total2, 0);

                    boolean[] roiMask = roi.MakeMask();

                    for (int y = 0, index = 0;
                         y < roi.dims.height;
                         ++y)
                    {
                        for (int x = 0;
                             x < roi.dims.width;
                             ++x, ++index)
                        {
                            if (roiMask[index])
                            {
                                int globalY = y + roi.origin.y;
                                int globalX = x + roi.origin.x;

                                if (globalY >= 0
                                    && globalY < dims.height
                                    && globalX >= 0
                                    && globalX < dims.width)
                                {
                                    int globalIndex = (globalY
                                                       * dims.width)
                                                      + globalX;

                                    if ((mask[globalIndex] & 0xff) > 0)
                                    {
                                        ++nPix;

                                        for (int chIndex = 0;
                                             chIndex < nChannels;
                                             ++chIndex)
                                        {
                                            double value = (double) frame.channelSelector.
                                                    get(CHANNEL_OUTPUT_ORDER[chIndex])[globalIndex];

                                            total[chIndex] += value;
                                            total2[chIndex] += (value * value);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    for (int chIndex = 0;
                         chIndex < nChannels;
                         ++chIndex)
                    {
                        if (nPix > 0)
                        {

                            Statistics stat = new Statistics();
                            stat.channelName = CHANNEL_OUTPUT_ORDER[chIndex];

                            stat.total = total[chIndex];
                            stat.mean = stat.total / nPix;
                            stat.sd = Math.sqrt((total2[chIndex] / nPix)
                                                - (stat.mean * stat.mean));

                            stats[chIndex][myTaskIndex] = stat;
                        }
                        else
                        {
                            stats[chIndex][myTaskIndex] = null;
                        }
                    }
                };

                tasks[taskIndex++] = multithreader.submit(task);
            }

            for (Future<?> future : tasks)
            {
                try
                {
                    future.get();
                }
                catch (InterruptedException | ExecutionException ex)
                {
                    System.err.println(ex.getMessage());
                }
            }
        }
    }

    private
            class ChannelLookup
    {

        HashMap<String, Integer> Channels;

        public
                ChannelLookup(MetadataStore metaStore)
        {
            OMEXMLMetadataImpl meta = (OMEXMLMetadataImpl) metaStore;
            Channels = new HashMap<>();

            for (int index = 0;
                 index < meta.getChannelCount(0);
                 ++index)
            {
                Channels.put(meta.getChannelName(0, index), index);
            }
        }

        public
                int GetChannelIndex(String name) throws FormatException
        {
            Integer val = Channels.get(name);

            if (val == null)
            {
                throw new FormatException("Channel name: " + name
                                          + " not found");
            }

            return val;
        }
    }

    /*
        ////////    movie output classes   START    ////////
     */
    interface MovieFrameUpdater
    {

        void FrameUpdate(BufferedImage img);

        void Close();
    }

    class MovieFrameGenerator implements ContourTracker.ContourUpdate
    {

        private
                MovieFrameUpdater updater;

        private
                Future<?> updaterJob;

        private final
                AtomicInteger nWriting;

        private final
                AffineTransformOp op;

        private final
                AffineTransform at;

        private final
                Dims inputDims;

        private final
                Dims outputDims;

        private
                BufferedImage backgroundImg;

        private final
                HashMap<Contour, GeneralPath> paths;

        private
                GeneralPath[] shapes;

        private
                int timepoint;

        private final
                String name;

        private final
                int detailedContourUpdateInterval;

        MovieFrameGenerator(byte[] pixels, Dims dims, double scaleFactor,
                            String name, int detailedContourUpdateInterval)
        {
            timepoint = -1;
            updater = null;
            this.name = name;
            updaterJob = null;
            this.detailedContourUpdateInterval = detailedContourUpdateInterval;
            nWriting = new AtomicInteger(0);

            at = AffineTransform.
                    getScaleInstance(scaleFactor, scaleFactor);
            op = new AffineTransformOp(at, AffineTransformOp.TYPE_BICUBIC);
            inputDims = dims;

            paths = new HashMap<>();
            shapes = new GeneralPath[0];

            ChangeImage(pixels);
            outputDims = new Dims(backgroundImg.getWidth(),
                                  backgroundImg.getHeight());
        }

        Dims GetOutputDims()
        {
            return outputDims;
        }

        void SetUpdater(MovieFrameUpdater updater)
        {
            this.updater = updater;
        }

        private
                void ChangeImage(byte[] pixels)
        {
            BufferedImage img = new BufferedImage(inputDims.width, inputDims.height,
                                                  BufferedImage.TYPE_BYTE_GRAY);
            System.arraycopy(pixels, 0, ((DataBufferByte) img.getRaster().
                                         getDataBuffer()).getData(), 0, pixels.length);
            backgroundImg = op.filter(img, null);

            ++timepoint;
            RemoveAllContours();
        }

        private
                GeneralPath GeneratePath(Contour contour)
        {

            GeneralPath shape = new GeneralPath();
            Vector2[] points = contour.points.toArray(new Vector2[0]);

            shape.moveTo(points[0].x, points[0].y);
            for (int index = 1;
                 index < points.length;
                 ++index)
            {
                shape.lineTo(points[index].x, points[index].y);
            }
            shape.closePath();

            shape.transform(at);

            return shape;
        }

        private
                void Update(boolean updateOverride)
        {
            if (updater != null)
            {
                boolean detailedUpdate = timepoint >= 0 && (timepoint
                                                            % detailedContourUpdateInterval) == 0;

                if (detailedUpdate != updateOverride)
                {
                    synchronized (nWriting)
                    {
                        nWriting.incrementAndGet();
                    }

                    SetPaths(paths.values().
                            toArray(new GeneralPath[0]));

                    BufferedImage outputImg = new BufferedImage(
                            outputDims.width,
                            outputDims.height,
                            BufferedImage.TYPE_3BYTE_BGR);

                    Paint(outputImg.getGraphics(), detailedUpdate);

                    if (updaterJob != null)
                    {
                        try
                        {
                            updaterJob.get();
                        }
                        catch (InterruptedException | ExecutionException ex)
                        {
                            Logger.getLogger(ActiveContourTracker.class.getName()).log(Level.SEVERE, null, ex);
                        }
                    }

                    updaterJob = multithreader.submit(()
                            ->
                    {
                        updater.FrameUpdate(outputImg);

                        synchronized (nWriting)
                        {
                            nWriting.decrementAndGet();
                        }
                    });
                }
            }
        }

        public
                void PushUpdate()
        {
            Update(true);
        }

        @Override
        public
                void Refresh()
        {
            Update(false);
        }

        @Override
        public
                void AddContour(Contour contour)
        {
            paths.put(contour, GeneratePath(contour));
        }

        @Override
        public
                void UpdateContour(Contour contour)
        {
            AddContour(contour);
        }

        @Override
        public
                void RemoveContour(Contour contour)
        {
            paths.remove(contour);
        }

        void RemoveAllContours()
        {
            paths.clear();
        }

        void SetPaths(GeneralPath[] shapes)
        {
            this.shapes = shapes;
        }

        GeneralPath[] GetContourShapes()
        {
            return shapes;
        }

        private
                void PaintBackground(java.awt.Graphics2D g2)
        {
            g2.drawImage(backgroundImg, new AffineTransform(), null);
        }

        private
                void DrawShapes(java.awt.Graphics2D g2)
        {
            for (GeneralPath shape : shapes)
            {
                g2.setPaint(java.awt.Color.RED);
                g2.draw(shape);
                g2.setPaint(java.awt.Color.GREEN);
                g2.fill(shape);
            }
        }

        public
                void Paint(java.awt.Graphics g, boolean detailedUpdate)
        {
            java.awt.Graphics2D g2 = (java.awt.Graphics2D) g;
            g2.setRenderingHint(java.awt.RenderingHints.KEY_ANTIALIASING,
                                java.awt.RenderingHints.VALUE_ANTIALIAS_ON);
            PaintBackground(g2);

            g2.setStroke(new java.awt.BasicStroke(2));
            g2.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 0.50f));
            DrawShapes(g2);

            String tag = name + "\nt: " + timepoint;
            if (detailedUpdate)
            {
                tag += ("\nContours: " + shapes.length);
            }

            g2.setComposite(java.awt.AlphaComposite.getInstance(java.awt.AlphaComposite.SRC_OVER, 1f));
            g2.setFont(new Font("Verdana", Font.BOLD, 16 * outputDims.width / 700));
            g2.setColor(java.awt.Color.WHITE);
            DrawString(g2, tag, outputDims.width / 10, outputDims.height / 10);
        }

        private
                void DrawString(java.awt.Graphics2D g2, String text, float x, float y)
        {
            for (String line : text.split("\n"))
            {
                g2.drawString(line, x, y += g2.getFontMetrics().getHeight());
            }
        }

        void Close() throws InterruptedException
        {
            if (updater != null)
            {
                int n;
                int timeOutMills = 600000;
                int t = 0;
                do
                {
                    synchronized (nWriting)
                    {
                        n = nWriting.intValue();
                    }
                    if (n > 0)
                    {
                        Thread.sleep(1);
                        ++t;
                    }
                }
                while (n > 0 && t < timeOutMills);

                if (n > 0)
                {
                    throw new InterruptedException("Timed out waiting for "
                                                   + "movie writer to close");
                }

                updater.Close();
            }
        }
    }

    private
            class SimpleMovieFrameGenerator
    {

        private
                class ContourUpdateStore implements ContourTracker.ContourUpdate
        {

            private final
                    HashMap<Contour, GeneralPath> contours;

            ContourUpdateStore()
            {
                contours = new HashMap<>();
            }

            private
                    GeneralPath GeneratePath(Contour contour)
            {

                GeneralPath shape = new GeneralPath();
                Vector2[] points = contour.points.toArray(new Vector2[0]);

                shape.moveTo(points[0].x, points[0].y);
                for (int index = 1;
                     index < points.length;
                     ++index)
                {
                    shape.lineTo(points[index].x, points[index].y);
                }
                shape.closePath();

                shape.transform(at);

                return shape;
            }

            void RemoveAllContours()
            {
                contours.clear();
            }

            @Override
            public
                    void AddContour(Contour contour)
            {
                contours.put(contour, GeneratePath(contour));
            }

            @Override
            public
                    void RemoveContour(Contour contour)
            {
                contours.remove(contour);
            }

            @Override
            public
                    void UpdateContour(Contour contour)
            {
                AddContour(contour);
            }

            @Override
            public
                    void Refresh()
            {
                // Do nothing
            }

            GeneralPath[] GetContourPaths()
            {
                return contours.values().toArray(new GeneralPath[0]);
            }
        }

        private
                MovieFrameUpdater updater;

        private
                Future<?> updaterJob;

        private final
                AtomicInteger nWriting;

        private final
                AffineTransformOp op;

        private final
                AffineTransform at;

        private final
                Dims inputDims;

        private final
                Dims pannelOutputDims;

        private final
                Dims movieOutputDims;

        private
                int timepoint;

        private final
                ContourUpdateStore nucleiContours;

        private final
                ContourUpdateStore wholeCellContours;

        SimpleMovieFrameGenerator(Dims dims, double scaleFactor)
        {
            timepoint = -1;
            inputDims = dims;
            updaterJob = null;
            nWriting = new AtomicInteger(0);

            at = AffineTransform.getScaleInstance(scaleFactor, scaleFactor);
            op = new AffineTransformOp(at, AffineTransformOp.TYPE_BICUBIC);

            BufferedImage testImg = op.filter(new BufferedImage(inputDims.width,
                                                                inputDims.height, BufferedImage.TYPE_BYTE_GRAY), null);

            pannelOutputDims = new Dims(testImg.getWidth(), testImg.getHeight());

            movieOutputDims = new Dims(3 * pannelOutputDims.width,
                                       2 * pannelOutputDims.height);

            nucleiContours = new ContourUpdateStore();
            wholeCellContours = new ContourUpdateStore();
        }

        void SetUpdater(MovieFrameUpdater updater)
        {
            this.updater = updater;
        }

        Dims GetMovieFrameDims()
        {
            return movieOutputDims;
        }

        ContourUpdateStore GetNucleiContours()
        {
            return nucleiContours;
        }

        ContourUpdateStore GetWholeCellContours()
        {
            return wholeCellContours;
        }

        void Close() throws InterruptedException
        {
            if (updater != null)
            {
                int n;
                int timeOutMills = 600000;
                int t = 0;
                do
                {
                    synchronized (nWriting)
                    {
                        n = nWriting.intValue();
                    }
                    if (n > 0)
                    {
                        Thread.sleep(1);
                        ++t;
                    }
                }
                while (n > 0 && t < timeOutMills);

                if (n > 0)
                {
                    throw new InterruptedException("Timed out waiting for "
                                                   + "movie writer to close");
                }

                updater.Close();
            }
        }

        private
                short[] NormalisePixels(short[] in)
        {
            float min = Float.MAX_VALUE;
            float max = 0;

            float shortMax = Short.MAX_VALUE & 0xffff;

            for (int index = 0;
                 index < in.length;
                 ++index)
            {
                float value = in[index] & 0xffff;
                if (value < min)
                {
                    min = value;
                }
                if (value > max)
                {
                    max = value;
                }
            }
            shortMax /= (max - min);

            short[] out = new short[in.length];
            for (int index = 0;
                 index < in.length;
                 ++index)
            {
                out[index] = (short) (((in[index] & 0xffff) - min) * shortMax);
            }

            return out;
        }

        void SetFrame(Frame frame)
        {
            if (updater != null)
            {

                synchronized (nWriting)
                {
                    nWriting.incrementAndGet();
                }

                ++timepoint;

                int nPix = inputDims.width * inputDims.height;

                BufferedImage[] backgroundImgs
                                = new BufferedImage[CHANNEL_OUTPUT_ORDER.length + 2];
                for (int index = 0;
                     index < backgroundImgs.length;
                     ++index)
                {
                    backgroundImgs[index] = new BufferedImage(pannelOutputDims.width,
                                                              pannelOutputDims.height, BufferedImage.TYPE_INT_ARGB);
                }

                BufferedImage img = new BufferedImage(inputDims.width,
                                                      inputDims.height, BufferedImage.TYPE_BYTE_GRAY);

                System.arraycopy(frame.fucciMask, 0, ((DataBufferByte) img.getRaster().
                                                      getDataBuffer()).getData(), 0, nPix);
                op.filter(img, backgroundImgs[0]);

                System.arraycopy(frame.redoxMask, 0, ((DataBufferByte) img.getRaster().
                                                      getDataBuffer()).getData(), 0, nPix);
                op.filter(img, backgroundImgs[1]);

                img = new BufferedImage(inputDims.width,
                                        inputDims.height, BufferedImage.TYPE_USHORT_GRAY);
                for (int index = 0;
                     index < CHANNEL_OUTPUT_ORDER.length;
                     ++index)
                {
                    System.arraycopy(NormalisePixels(frame.channelSelector.get(CHANNEL_OUTPUT_ORDER[index])),
                                     0, ((DataBufferUShort) img.getRaster().getDataBuffer()).getData(), 0, nPix);
                    op.filter(img, backgroundImgs[index + 2]);
                }

                GeneralPath[] wholeCellPaths = wholeCellContours.GetContourPaths();
                wholeCellContours.RemoveAllContours();
                GeneralPath[] nucleiPaths = nucleiContours.GetContourPaths();
                nucleiContours.RemoveAllContours();

                Color[] wholeCellColours = new Color[]
                {
                    Color.YELLOW, Color.BLUE
                };
                Color[] nucleiColours = new Color[]
                {
                    Color.RED, Color.GREEN
                };

                float x = pannelOutputDims.width / 10;
                float y = pannelOutputDims.height / 10;

                for (int index = 0;
                     index < backgroundImgs.length;
                     ++index)
                {
                    Graphics2D g2 = (Graphics2D) backgroundImgs[index].getGraphics();

                    if (index != 0)
                    {
                        DrawShapes(g2, wholeCellPaths,
                                   wholeCellColours[0], wholeCellColours[1]);
                    }

                    if (index != 1)
                    {
                        DrawShapes(g2, nucleiPaths,
                                   nucleiColours[0], nucleiColours[1]);
                    }

                    switch (index)
                    {
                        case 0:
                            y = DrawString(g2, "t: " + timepoint, x, y, Color.WHITE);
                            DrawString(g2, "Nuclei Mask", x, y, nucleiColours[1]);
                            break;
                        case 1:
                            DrawString(g2, "Whole Cell Mask", x, y, wholeCellColours[1]);
                            break;
                        default:
                            DrawString(g2, CHANNEL_OUTPUT_ORDER[index - 2], x, y, Color.WHITE);
                            break;
                    }

                    g2.setStroke(new BasicStroke(3));
                    g2.setColor(Color.ORANGE);
                    g2.drawRect(0, 0, pannelOutputDims.width, pannelOutputDims.height);
                    g2.dispose();
                }

                BufferedImage movieFrame = new BufferedImage(movieOutputDims.width,
                                                             movieOutputDims.height, BufferedImage.TYPE_3BYTE_BGR);

                Graphics2D g2 = (Graphics2D) movieFrame.getGraphics();

                g2.drawImage(backgroundImgs[0], 0, 0, null);
                g2.drawImage(backgroundImgs[1], 0, pannelOutputDims.height, null);
                g2.drawImage(backgroundImgs[2], pannelOutputDims.width, 0, null);
                g2.drawImage(backgroundImgs[3], pannelOutputDims.width,
                             pannelOutputDims.height, null);
                g2.drawImage(backgroundImgs[4], 2 * pannelOutputDims.width, 0, null);
                g2.drawImage(backgroundImgs[5], 2 * pannelOutputDims.width,
                             pannelOutputDims.height, null);

                g2.dispose();

                if (updaterJob != null)
                {
                    try
                    {
                        updaterJob.get();
                    }
                    catch (InterruptedException | ExecutionException ex)
                    {
                        Logger.getLogger(ActiveContourTracker.class.getName()).
                                log(Level.SEVERE, null, ex);
                    }
                }

                updaterJob = multithreader.submit(()
                        ->
                {
                    updater.FrameUpdate(movieFrame);

                    synchronized (nWriting)
                    {
                        nWriting.decrementAndGet();
                    }
                });
            }
        }

        private
                float DrawString(Graphics2D g2, String text, float x, float y,
                                 Color c)
        {
            g2.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 1f));
            g2.setFont(new Font("Verdana", Font.BOLD, 16 * pannelOutputDims.width / 700));
            g2.setColor(c);
            for (String line : text.split("\n"))
            {
                g2.drawString(line, x, y += g2.getFontMetrics().getHeight());
            }

            return y;
        }

        private
                void DrawShapes(Graphics2D g2, GeneralPath[] shapes,
                                Color outline, Color fill)
        {
            g2.setStroke(new BasicStroke(2));
            g2.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER,
                                                       0.50f));
            for (GeneralPath shape : shapes)
            {
                g2.setPaint(outline);
                g2.draw(shape);
                g2.setPaint(fill);
                g2.fill(shape);
            }
        }

    }
    /*
        ////////    movie output classes    END     ////////
     */
}
