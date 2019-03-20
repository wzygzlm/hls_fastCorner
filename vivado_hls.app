<project xmlns="com.autoesl.autopilot.project" name="hls_fastCorner" top="fastCornerHW">
    <includePaths/>
    <libraryPaths/>
    <Simulation>
        <SimFlow name="csim"/>
    </Simulation>
    <files xmlns="">
        <file name="../src/sortTB.cpp" sc="0" tb="1" cflags=" "/>
        <file name="hls_fastCorner/src/sortHW.cpp" sc="0" tb="false" cflags=""/>
        <file name="hls_fastCorner/src/fast_detector.cpp" sc="0" tb="false" cflags=""/>
        <file name="hls_fastCorner/src/fastCorner.cpp" sc="0" tb="false" cflags=""/>
    </files>
    <solutions xmlns="">
        <solution name="solution1" status="active"/>
    </solutions>
</project>

