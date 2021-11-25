const path = require("path")
module.exports = {
    entry: './src/index.js',
    output:{
        path: path.resolve(__dirname, './wgscovplot/tmpl/js'),
        filename: "wgscovplot.bundle.js",
        library: "wgscovplot",
        libraryTarget: "umd",
        globalObject: "this"
    },
    module: {
        rules:[{
            test:/\.js$/,
            exclude:[/(node_modules)/,/(wgscovplot)/],
            use: "babel-loader"
        }],
    },
};