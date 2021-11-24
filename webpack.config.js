const path = require("path")
module.exports = {
    mode: "production",
    entry: './src/index.js',
    output:{
        path: path.resolve(__dirname, './wgscovplot/tmpl/js'),
        filename: "wgscovplot.js",
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