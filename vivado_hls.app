<project xmlns="com.autoesl.autopilot.project" name="hls_fastCorner" top="fastCornerInnerHW">
    <includePaths/>
    <libraryPaths/>
    <Simulation>
        <SimFlow name="csim" csimMode="0" lastCsimMode="0"/>
    </Simulation>
    <files xmlns="">
        <file name="../src/sortTB.cpp" sc="0" tb="1" cflags=" "/>
        <file name="hls_fastCorner/src/fast_detector.cpp" sc="0" tb="false" cflags=""/>
        <file name="hls_fastCorner/src/fastCorner.cpp" sc="0" tb="false" cflags=""/>
    </files>
    <solutions xmlns="">
        <solution name="solution1" status="active"/>
        <solution name="solution2" status="inactive"/>
    </solutions>
</project>

