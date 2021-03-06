/*
    A lot of this is copied from Humble-Video screenshot output demo
 */
package activecontourtracker;

/**
 *
 * @author edwardharry
 */
import java.awt.AWTException;
import java.awt.image.BufferedImage;
import java.io.IOException;

import io.humble.video.Codec;
import io.humble.video.Encoder;
import io.humble.video.MediaPacket;
import io.humble.video.MediaPicture;
import io.humble.video.Muxer;
import io.humble.video.MuxerFormat;
import io.humble.video.PixelFormat;
import io.humble.video.Rational;
import io.humble.video.awt.MediaPictureConverter;
import io.humble.video.awt.MediaPictureConverterFactory;

public
        class MovieCreator
{

    static
            Movie VP9Movie(String filename, int snapsPerSecond, Dims dims)
            throws AWTException, InterruptedException, IOException
    {
        return new Movie(filename, null, Codec.ID.CODEC_ID_VP9,
                         snapsPerSecond, dims);
    }

    public static
            class Movie implements ActiveContourTracker.MovieFrameUpdater
    {

        private final
                Encoder encoder;
        private final
                MediaPicture picture;
        private final
                MediaPacket packet;
        private final
                Muxer muxer;

        private
                int timepoint;

        Movie(String filename, String formatname,
              String codecname, int snapsPerSecond, Dims dims)
                throws AWTException, InterruptedException, IOException
        {
            timepoint = 0;

            final
                    Rational framerate = Rational.make(1, snapsPerSecond);

            /**
             * First we create a muxer using the passed in filename and
             * formatname if given.
             */
            muxer = Muxer.make(filename, null, formatname);

            /**
             * Now, we need to decide what type of codec to use to encode video.
             * Muxers have limited sets of codecs they can use. We're going to
             * pick the first one that works, or if the user supplied a codec
             * name, we're going to force-fit that in instead.
             */
            final
                    MuxerFormat format = muxer.getFormat();
            final
                    Codec codec;

            if (codecname != null)
            {
                codec = Codec.findEncodingCodecByName(codecname);
            }
            else
            {
                codec = Codec.findEncodingCodec(format.getDefaultVideoCodecId());
            }

            /**
             * Now that we know what codec, we need to create an encoder
             */
            encoder = Encoder.make(codec);

            /**
             * Video encoders need to know at a minimum: width height pixel
             * format Some also need to know frame-rate (older codecs that had a
             * fixed rate at which video files could be written needed this).
             * There are many other options you can set on an encoder, but we're
             * going to keep it simpler here.
             */
            encoder.setWidth(dims.width);
            encoder.setHeight(dims.height);
            // We are going to use 420P as the format because that's what most video formats these days use
            final
                    PixelFormat.Type pixelformat = PixelFormat.Type.PIX_FMT_YUV420P;
            encoder.setPixelFormat(pixelformat);
            encoder.setTimeBase(framerate);

            /**
             * An annoynace of some formats is that they need global (rather
             * than per-stream) headers, and in that case you have to tell the
             * encoder. And since Encoders are decoupled from Muxers, there is
             * no easy way to know this beyond
             */
            if (format.getFlag(MuxerFormat.Flag.GLOBAL_HEADER))
            {
                encoder.setFlag(Encoder.Flag.FLAG_GLOBAL_HEADER, true);
            }

            /**
             * Open the encoder.
             */
            encoder.open(null, null);

            /**
             * Add this stream to the muxer.
             */
            muxer.addNewStream(encoder);

            /**
             * And open the muxer for business.
             */
            muxer.open(null, null);

            /**
             * Next, we need to make sure we have the right MediaPicture format
             * objects to encode data with. Java (and most on-screen graphics
             * programs) use some variant of Red-Green-Blue image encoding
             * (a.k.a. RGB or BGR). Most video codecs use some variant of YCrCb
             * formatting. So we're going to have to convert. To do that, we'll
             * introduce a MediaPictureConverter object later. object.
             */
            picture = MediaPicture.make(
                    encoder.getWidth(),
                    encoder.getHeight(),
                    pixelformat);
            picture.setTimeBase(framerate);

            packet = MediaPacket.make();
        }

        Movie(String filename, String formatname,
              Codec.ID id, int snapsPerSecond, Dims dims)
                throws AWTException, InterruptedException, IOException
        {
            timepoint = 0;

            final
                    Rational framerate = Rational.make(1, snapsPerSecond);

            /**
             * First we create a muxer using the passed in filename and
             * formatname if given.
             */
            muxer = Muxer.make(filename, null, formatname);

            /**
             * Now, we need to decide what type of codec to use to encode video.
             * Muxers have limited sets of codecs they can use. We're going to
             * pick the first one that works, or if the user supplied a codec
             * name, we're going to force-fit that in instead.
             */
            final
                    MuxerFormat format = muxer.getFormat();
            final
                    Codec codec = Codec.findEncodingCodec(id);

            /**
             * Now that we know what codec, we need to create an encoder
             */
            encoder = Encoder.make(codec);

            /**
             * Video encoders need to know at a minimum: width height pixel
             * format Some also need to know frame-rate (older codecs that had a
             * fixed rate at which video files could be written needed this).
             * There are many other options you can set on an encoder, but we're
             * going to keep it simpler here.
             */
            encoder.setWidth(dims.width);
            encoder.setHeight(dims.height);
            // We are going to use 420P as the format because that's what most video formats these days use
            final
                    PixelFormat.Type pixelformat = PixelFormat.Type.PIX_FMT_YUV420P;
            encoder.setPixelFormat(pixelformat);
            encoder.setTimeBase(framerate);

            /**
             * An annoynace of some formats is that they need global (rather
             * than per-stream) headers, and in that case you have to tell the
             * encoder. And since Encoders are decoupled from Muxers, there is
             * no easy way to know this beyond
             */
            if (format.getFlag(MuxerFormat.Flag.GLOBAL_HEADER))
            {
                encoder.setFlag(Encoder.Flag.FLAG_GLOBAL_HEADER, true);
            }

            /**
             * Open the encoder.
             */
            encoder.open(null, null);

            /**
             * Add this stream to the muxer.
             */
            muxer.addNewStream(encoder);

            /**
             * And open the muxer for business.
             */
            muxer.open(null, null);

            /**
             * Next, we need to make sure we have the right MediaPicture format
             * objects to encode data with. Java (and most on-screen graphics
             * programs) use some variant of Red-Green-Blue image encoding
             * (a.k.a. RGB or BGR). Most video codecs use some variant of YCrCb
             * formatting. So we're going to have to convert. To do that, we'll
             * introduce a MediaPictureConverter object later. object.
             */
            picture = MediaPicture.make(
                    encoder.getWidth(),
                    encoder.getHeight(),
                    pixelformat);
            picture.setTimeBase(framerate);

            packet = MediaPacket.make();
        }

        private
                void WriteFrame(BufferedImage imgIn)
        {
            /**
             * Make the screen capture && convert image to TYPE_3BYTE_BGR
             */
            final
                    BufferedImage imgRec = ConvertToType(imgIn,
                                                         BufferedImage.TYPE_3BYTE_BGR);

            /**
             * This is LIKELY not in YUV420P format, so we're going to convert
             * it using some handy utilities.
             */
            MediaPictureConverter converter = MediaPictureConverterFactory.
                    createConverter(imgRec, picture);
            converter.toPicture(picture, imgRec, timepoint++);

            synchronized (muxer)
            {
                do
                {
                    encoder.encode(packet, picture);
                    if (packet.isComplete())
                    {
                        muxer.write(packet, false);
                    }
                }
                while (packet.isComplete());
            }

        }

        @Override
        public
                void Close()
        {
            synchronized (muxer)
            {
                /**
                 * Encoders, like decoders, sometimes cache pictures so it can
                 * do the right key-frame optimizations. So, they need to be
                 * flushed as well. As with the decoders, the convention is to
                 * pass in a null input until the output is not complete.
                 */
                do
                {
                    encoder.encode(packet, null);
                    if (packet.isComplete())
                    {
                        muxer.write(packet, false);
                    }
                }
                while (packet.isComplete());

                /**
                 * Finally, let's clean up after ourselves.
                 */
                muxer.close();
            }
        }

        @Override
        public
                void FrameUpdate(BufferedImage img)
        {
            WriteFrame(img);
        }
    }

    /**
     * Convert a {@link BufferedImage} of any type, to {@link BufferedImage} of
     * a specified type. If the source image is the same type as the target
     * type, then original image is returned, otherwise new image of the correct
     * type is created and the content of the source image is copied into the
     * new image.
     *
     * @param sourceImage the image to be converted
     * @param targetType the desired BufferedImage type
     *
     * @return a BufferedImage of the specified target type.
     *
     * @see BufferedImage
     */
    public static
            BufferedImage ConvertToType(BufferedImage sourceImage,
                                        int targetType)
    {
        BufferedImage image;

        // if the source image is already the target type, return the source image
        if (sourceImage.getType() == targetType)
        {
            image = sourceImage;
        }

        // otherwise create a new image of the target type and draw the new
        // image
        else
        {
            image = new BufferedImage(sourceImage.getWidth(),
                                      sourceImage.getHeight(), targetType);
            image.getGraphics().drawImage(sourceImage, 0, 0, null);
        }

        return image;
    }

}
